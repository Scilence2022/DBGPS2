#!/usr/bin/env python3
"""Explicit design-assisted supplementation controls; never sequencing-only recovery."""
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time
import zipfile

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'input/antkowiak2020'
OUT = DATA / 'results/file3_supplement_simulation'
PREVIOUS = DATA / 'results/block_decode'
SMALL = DATA / 'results/block12_small_k'
B = 16383
START = 11 * B
N = 12 * B
BASE4 = str.maketrans('ACGT', '0123')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    rows = {}
    with path.open() as f:
        for line in f:
            field, dna = line.rstrip().split('\t')
            i = int(field)
            assert i not in rows and 0 <= i < N and len(dna) == 60 and set(dna) <= set('ACGT')
            assert (int(dna[42:54].translate(BASE4), 4) ^ 1071068) == i
            rows[i] = dna
    return rows


def save(path, rows):
    with path.open('w') as f:
        for i, dna in sorted(rows.items()):
            f.write(f'{i}\t{dna}\n')


def score(rows, truth, block):
    assert all((block - 1) * B <= i < block * B for i in rows)
    correct = 0
    errors = [0] * 6
    for i, dna in rows.items():
        if dna == truth[i]:
            correct += 1
        else:
            delta = int(dna[:42].translate(BASE4), 4) ^ int(truth[i][:42].translate(BASE4), 4)
            assert delta
            for lane in range(6):
                errors[lane] += bool((delta >> (14 * lane)) & 16383)
    erasures = B - len(rows)
    budgets = [erasures + 2 * e for e in errors]
    return dict(block=block, emitted=len(rows), correct=correct, wrong=len(rows)-correct,
                erasures=erasures, lane_errors=errors, lane_budgets=budgets,
                max_lane_budget=max(budgets), within_rs_bound=max(budgets) <= 5406)


def decode(case, reference):
    folder = Path(case['folder'])
    output = folder / 'simulated_recovered.zip'
    command = [str(DATA / 'bin/rs_driver_fast'), 'decode', case['combined_input'], str(output), '12']
    print('START DECODE', case['name'], flush=True)
    start = time.monotonic()
    with (folder / 'outer_decode.log').open('w') as log:
        process = subprocess.run(command, stdout=log, stderr=log, timeout=1200)
    result = dict(command=command, exit_code=process.returncode, seconds=time.monotonic()-start,
                  path=str(output), simulated=True)
    if output.exists():
        result.update(bytes=output.stat().st_size, sha256=sha(output), byte_identical_to_reference=output.read_bytes() == reference.read_bytes())
        try:
            with zipfile.ZipFile(output) as archive:
                result['zip_crc_error'] = archive.testzip()
                result['zip_entries'] = [dict(name=i.filename, bytes=i.file_size, crc32=f'{i.CRC:08x}') for i in archive.infolist()]
        except (zipfile.BadZipFile, RuntimeError, EOFError) as error:
            result['zip_error'] = str(error)
    case['decode'] = result
    (folder / 'result.json').write_text(json.dumps(case, indent=2))
    print('FINISHED', case['name'], json.dumps({key: result.get(key) for key in
          ['exit_code', 'seconds', 'bytes', 'sha256', 'byte_identical_to_reference', 'zip_crc_error']}), flush=True)
    assert process.returncode == 0 and result.get('byte_identical_to_reference') and result.get('zip_crc_error') is None
    assert 'zip_error' not in result
    return case


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    truth_path = DATA / 'sources/File3_synthesized.fasta'
    truth = [s for s in truth_path.read_text().splitlines() if s and not s.startswith('>')]
    assert len(truth) == N
    reference = DATA / 'results/reference_file3.zip'
    assert sha(reference) == 'a7d53ade64bab3f60252e75b28e1741d85fdd7b5b393ffddc95fd141466d8cbb'
    previous = json.loads((PREVIOUS / 'summary.json').read_text())
    blocks = []
    prefix = bytearray()
    for b, chosen in enumerate(previous['chosen_per_block'][:11], 1):
        path = Path(chosen['selected_strands'])
        rows = load(path)
        measured = score(rows, truth, b)
        assert measured['within_rs_bound']
        assert measured['correct'] == chosen['correct'] and measured['wrong'] == chosen['errors']
        assert measured['lane_budgets'] == chosen['lane_budgets']
        blocks.append(dict(**measured, input=str(path), sha256=sha(path), changed=False))
        prefix.extend(path.read_bytes())
    strict_path = SMALL / 'best_budget.strands'
    strict = load(strict_path)
    before = score(strict, truth, 12)
    assert (before['correct'], before['wrong'], before['erasures']) == (9973, 15, 6395)
    with (SMALL / 'per_index.tsv').open() as f:
        indices = [int(row['index']) for row in csv.DictReader(f, delimiter='\t') if int(row['combined_success'])]
    assert len(set(indices)) == 10659 and all(START <= i < N for i in indices)
    oracle = {i: truth[i] for i in indices}
    cases = []
    for name, base, count, description in [
        ('actual_plus_1019', strict, 1019, 'Actual k16/c1/r3.5 strict output, retaining all 15 wrong strands; inject correct design strands only at missing indices.'),
        ('oracle_union_plus_318', oracle, 318, 'Truth-select the 10659 indices with an exact candidate under at least one of 353 parameters, then inject correct design strands at 318 missing indices.'),
    ]:
        folder = OUT / name
        folder.mkdir(exist_ok=True)
        missing = [i for i in range(START, N) if i not in base]
        added = missing[:count]
        supplemented = dict(base)
        supplemented.update({i: truth[i] for i in added})
        assert len(added) == count and all(supplemented[i] == dna for i, dna in base.items())
        assert all(i not in base for i in added)
        original = score(base, truth, 12)
        after = score(supplemented, truth, 12)
        assert after['correct'] == original['correct'] + count and after['wrong'] == original['wrong']
        assert after['lane_errors'] == original['lane_errors'] and after['max_lane_budget'] == 5406
        # One fewer supplement leaves the worst lane one unit beyond the bound.
        one_fewer = dict(supplemented)
        del one_fewer[added[-1]]
        fewer_score = score(one_fewer, truth, 12)
        assert fewer_score['max_lane_budget'] == 5407
        save(folder / 'block12.before.strands', base)
        save(folder / 'block12.supplemented.strands', supplemented)
        with (folder / 'injected_design_strands.tsv').open('w') as f:
            f.write('index\tdna\tprovenance\n')
            for i in added:
                f.write(f'{i}\t{truth[i]}\tsimulated_design_injection\n')
        combined = folder / 'all12blocks.strands'
        combined.write_bytes(prefix + (folder / 'block12.supplemented.strands').read_bytes())
        parsed = load(combined)
        assert len(parsed) == sum(b['emitted'] for b in blocks) + len(supplemented)
        assert combined.read_bytes()[:len(prefix)] == prefix
        cases.append(dict(name=name, folder=str(folder), interpretation=description, simulated=True,
                          before=original, injected_correct=count, injection_rule='Ascending missing numeric index; never replace an existing strand.',
                          injected_indices=added, after=after, one_fewer_would_have=fewer_score,
                          combined_input=str(combined), combined_input_sha256=sha(combined),
                          initial_wrong_strands_preserved=True, first_11_blocks_unchanged=True))
    provenance = dict(simulated=True, reference_assisted_read_partition=True, source_truth=str(truth_path), source_truth_sha256=sha(truth_path),
                      decoder=str(DATA / 'bin/rs_driver_fast'), decoder_sha256=sha(DATA / 'bin/rs_driver_fast'),
                      reference_archive=str(reference), reference_archive_sha256=sha(reference),
                      first_11_blocks=blocks, cases=cases,
                      limitation='Design-injected and oracle-selected strands are not recovered from sequencing. Success demonstrates the outer decoder under the stated hypothetical input.')
    (OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        results = list(pool.map(lambda case: decode(case, reference), cases))
    provenance['cases'] = results
    provenance['all_cases_byte_identical'] = all(case['decode']['byte_identical_to_reference'] for case in results)
    (OUT / 'summary.json').write_text(json.dumps(provenance, indent=2))
    lines = ['# File 3 补入正确片段后的完整解码模拟', '',
             '**这是使用设计真值补片的模拟实验，不代表仅凭测序数据已恢复 File 3。**', '',
             '前 11 块保留已有实际唯一输出，逐字节不变；第 12 块仅在缺失 index 处补入完整正确设计片段，不替换任何已有片段。补片按缺失 index 升序选择。', '',
             '| 情形 | 补入正确片段 | 第12块正确 | 第12块错误 | 擦除 | 最差 2E+S | 完整文件一致 |',
             '|---|---:|---:|---:|---:|---:|---|']
    for case in results:
        a = case['after']
        lines.append(f"| {case['name']} | {case['injected_correct']} | {a['correct']} | {a['wrong']} | {a['erasures']} | {a['max_lane_budget']} | {case['decode']['byte_identical_to_reference']} |")
    lines += ['', 'actual_plus_1019：基于 k=16、c=1、r=3.5 的 strict 实际输出，原有 9973 条正确和 15 条错误全部保留。补后 10992 条正确 + 15 条错误，共 11007 条；6 路错误数不变，外码预算分别为 ' + str(results[0]['after']['lane_budgets']) + '。', '',
              'oracle_union_plus_318：先按真值选出 353 组参数中至少一次存在完整正确候选的 10659 个 index；该候选并集需要真值挑选。再补入 318 条缺失片段，得到 10977 条正确、0 条错误。', '',
              '两种情形都用原 RS 外码解码器 rs_driver_fast 对全部 12 块完整解码。解码器仅接收补片后的序列文件和块数，不读取设计序列或 ZIP 校验基准。', '',
              f"恢复 ZIP 均为 {reference.stat().st_size} 字节，与已有正确设计解码基准逐字节一致；所有 ZIP 条目通过 CRC 校验。SHA-256：`{sha(reference)}`。", '',
              '每种情形少补 1 条时，最差预算为 5407，超过上限 5406；这里只做预算计算，未对这些超界输入运行完整解码。', '',
              '补片清单见各情形目录的 injected_design_strands.tsv；完整输入见 all12blocks.strands；解码日志见 outer_decode.log；恢复文件见 simulated_recovered.zip；详细来源及核验结果见 summary.json。', '']
    (OUT / 'REPORT.zh-CN.md').write_text('\n'.join(lines))
    files = [OUT / 'summary.json', OUT / 'provenance.json', OUT / 'REPORT.zh-CN.md']
    for case in results:
        files.extend(p for p in Path(case['folder']).iterdir() if p.is_file())
    (OUT / 'validation.json').write_text(json.dumps(dict(source_sha256=sha(Path(__file__)), all_cases_byte_identical=True,
        artifact_sha256={str(path.relative_to(OUT)): sha(path) for path in files}), indent=2))
    print('SUPPLEMENT SIMULATION COMPLETE: BOTH FULL ARCHIVES EXACT', flush=True)


if __name__ == '__main__':
    main()
