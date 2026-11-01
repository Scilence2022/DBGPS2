# 解码代码包说明

这是论文 nbt.4079 的 DNA 片段恢复实验代码包，全部统计均为外层 RS 纠错前结果。
完整命令及文件说明见 README.md。建议使用 Python 3.10 及以上、C 编译器、zlib 和 pthreads。

请解压到不含空白字符的路径；历史C交互协议以空白分隔文件路径。

## 运行

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_suite.py --stage check
python run_suite.py --stage smoke
python validate_package.py
python run_suite.py --stage all
```

smoke 只运行包内20条片段的无误和替换错误对照；all 下载固定版本数据后运行全部已完成阶段。
允许使用 --stage 分阶段执行，顺序见英文说明。复现会使用解压目录，不需要原工作树或Git仓库。
代码包未包含原始大数据、完整参考集、逐候选大TSV或恢复FASTA，下载脚本会获取原始数据。
建议预留至少15 GB工作磁盘；此为预留建议，不是实测峰值。

historical_scripts 保留原始脚本；可运行副本只将两个脚本的Git版本查询替换为源码快照校验。
算法、参数和评分逻辑未改动。evidence 为原始小型记录，其中保留原机器绝对路径；
不要将这些记录复制到运行结果目录冒充可续跑状态，完整重跑会生成新路径。

480个初筛配置使用全量reads建图，只有查询种子按1/128采样；随后9个配置全量确认，
再运行19个正向残余种子配置。早期JSON按14 nt前缀分组，最终统计按13 nt编码头分组。
全局格式标定和调参使用参考序列，属于同数据集选优，未完成独立留出验证。
代码包保留原LICENSE；没有提交任何手稿或代码变更。
