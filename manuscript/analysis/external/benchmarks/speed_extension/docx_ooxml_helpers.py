"""Patch the integrated manuscript without rebuilding its existing OOXML."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

from lxml import etree as E
from PIL import Image

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[2]/'docx/DBGPS2_manuscript_with_figures_Gimpel2024.docx'
DESTINATION=SOURCE.with_name('DBGPS2_manuscript_speed_comparison.docx')
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'pic':'http://schemas.openxmlformats.org/drawingml/2006/picture',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W='{'+NS['w']+'}'


def element(tag,**attrs):
    return E.Element(W+tag,{W+k:str(v) for k,v in attrs.items()})


def text(node):
    return ''.join(node.itertext())


def replace_text(paragraph,value):
    for child in list(paragraph):
        if child.tag!=W+'pPr': paragraph.remove(child)
    run=E.SubElement(paragraph,W+'r');t=E.SubElement(run,W+'t')
    t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=value


def append_text(paragraph,value):
    run=E.SubElement(paragraph,W+'r');t=E.SubElement(run,W+'t')
    t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=value


def replace_run_text(paragraph,old,new):
    matches=[t for t in paragraph.findall('.//w:t',NS) if old in (t.text or '')]
    assert len(matches)==1,(old,len(matches))
    matches[0].text=matches[0].text.replace(old,new)


def paragraph(value,template):
    result=deepcopy(template);replace_text(result,value)
    first_run=template.find('w:r',NS)
    properties=first_run.find('w:rPr',NS) if first_run is not None else None
    if properties is not None:result.find('w:r',NS).insert(0,deepcopy(properties))
    return result


def table(rows):
    widths=[900,2900,1900,1400,900,1360]
    result=element('tbl');properties=E.SubElement(result,W+'tblPr')
    properties.extend([element('tblW',w=sum(widths),type='dxa'),element('jc',val='center'),element('tblLayout',type='fixed')])
    borders=E.SubElement(properties,W+'tblBorders')
    for side in ['top','left','bottom','right','insideH','insideV']:
        borders.append(element(side,val='single',sz=4,color='D9D9D9'))
    margins=E.SubElement(properties,W+'tblCellMar')
    for side in ['top','bottom']:margins.append(element(side,w=100,type='dxa'))
    for side in ['left','right']:margins.append(element(side,w=100,type='dxa'))
    grid=E.SubElement(result,W+'tblGrid')
    for width in widths:grid.append(element('gridCol',w=width))
    for row_index,values in enumerate(rows):
        row=E.SubElement(result,W+'tr');row_properties=E.SubElement(row,W+'trPr')
        row_properties.append(element('cantSplit'))
        if row_index==0:row_properties.append(element('tblHeader'))
        for column,(value,width) in enumerate(zip(values,widths)):
            cell=E.SubElement(row,W+'tc');cell_properties=E.SubElement(cell,W+'tcPr')
            cell_properties.extend([element('tcW',w=width,type='dxa'),element('vAlign',val='center')])
            if row_index==0:cell_properties.append(element('shd',fill='DCE6F1',val='clear'))
            elif row_index%2==0:cell_properties.append(element('shd',fill='F5F7F9',val='clear'))
            p=E.SubElement(cell,W+'p');ppr=E.SubElement(p,W+'pPr')
            ppr.extend([element('spacing',before=0,after=0,line=240,lineRule='auto'),element('jc',val='left' if column==1 else 'center')])
            run=E.SubElement(p,W+'r');rpr=E.SubElement(run,W+'rPr')
            rpr.extend([element('rFonts',ascii='Times New Roman',hAnsi='Times New Roman'),element('sz',val=20),element('color',val='000000')])
            if row_index==0:rpr.append(element('b'))
            t=E.SubElement(run,W+'t');t.text=str(value)
    return result

