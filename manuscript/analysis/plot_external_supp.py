"""Regenerate original manuscript figures from data/figure_data.json.
Requires Python 3.12, matplotlib 3.10.7 and numpy. All values are existing results.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

ROOT=Path(__file__).resolve().parents[1]
D=json.loads((ROOT/'data/source_data/figure_data.json').read_text())
OUT=ROOT/'figures'
BLUE='#2573B5';ORANGE='#D46A36';GREEN='#268478';RED='#B44D43';GRAY='#D5DADF';DARK='#26363E';PURPLE='#8272B2'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.2,'axes.titlesize':9.2,'axes.labelsize':8.3,'xtick.labelsize':7.6,'ytick.labelsize':7.6,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.65,'lines.linewidth':1.2,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})

def title(ax,letter,text):
 ax.set_title(text,loc='left',pad=10,fontweight='medium');ax.text(-0.13,1.055,letter,transform=ax.transAxes,fontweight='bold',fontsize=11,va='bottom')
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=400,bbox_inches='tight',pad_inches=.08)
 plt.close(fig)
def clean(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3,width=.6)
def diagram(ax):ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
def strand(ax,x,y,total,segments,width=.88,height=.13):
 cur=x
 for n,label,color in segments:
  w=width*n/total
  ax.add_patch(Rectangle((cur,y),w,height,facecolor=color,edgecolor='white',linewidth=.8))
  ax.text(cur+w/2,y+height/2,label,ha='center',va='center',fontsize=8,color='white' if color in [BLUE,GREEN,PURPLE,RED] else DARK)
  cur+=w

def fountain_bars(ax,letter='a'):
 title(ax,letter,'DNA Fountain code-valid outputs')
 f=D['fountain'];x=np.array([0,1]);w=.32
 for shift,key,color,label in [(-w/2,'merge',ORANGE,'Read merging'),(w/2,'DBGPS2',BLUE,'DBGPS2')]:
  vals=[f['full'][key],f['one_tile'][key]];ax.bar(x+shift,vals,w,color=color,label=label)
  for xx,v in zip(x+shift,vals):ax.text(xx,v+1000,f'{v:,}',ha='center',fontsize=7.7,rotation=0)
 ax.set_xticks(x,['Full depth','One tile']);ax.set_ylabel('Unique RS-valid strands');ax.set_ylim(0,95000);ax.set_yticks([0,20000,40000,60000,80000])
 ax.yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'))
 ax.legend(frameon=False,fontsize=7.5,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.02),handlelength=1)
 clean(ax)

fig=plt.figure(figsize=(7.35,5.2));gs=fig.add_gridspec(2,2,hspace=.6,wspace=.53)
fountain_bars(fig.add_subplot(gs[0,0]))
ax=fig.add_subplot(gs[0,1]);title(ax,'b','Photolithographic outputs')
x=np.arange(3);correct=np.array([v['correct']/v['designed']*100 for v in D['photolithographic']]);wrong=np.array([v['wrong']/v['designed']*100 for v in D['photolithographic']]);missing=100-correct-wrong
ax.bar(x,correct,.55,color=BLUE,label='Exact');ax.bar(x,wrong,.55,bottom=correct,color=RED,label='Wrong');ax.bar(x,missing,.55,bottom=correct+wrong,color=GRAY,label='Missing')
for j,v in enumerate(D['photolithographic']):
 ax.text(j,correct[j]/2,f"{v['recovery_pct']:.1f}%",ha='center',color='white',fontweight='bold',fontsize=8)
 ax.text(j,104,f"{v['wrong']} wrong",ha='center',fontsize=7.4)
ax.set_xticks(x,['File 1','File 2','File 3']);ax.set_ylabel('Designed indices (%)');ax.set_ylim(0,130);ax.set_yticks([0,25,50,75,100]);ax.legend(frameon=False,ncol=3,fontsize=7.2,loc='upper center',bbox_to_anchor=(.5,1.05),handlelength=.8,columnspacing=.8);clean(ax)
ax=fig.add_subplot(gs[1,0]);title(ax,'c','Organick recovery and errors')
colors=[GRAY,DARK,ORANGE,BLUE,GREEN,PURPLE];offsets=[(-10,6),(6,-13),(-13,-2),(-12,-1),(6,-2),(6,5)]
for row,c,off in zip(D['organick'],colors,offsets):
 ax.scatter(row['recovery_pct'],row['wrong'],s=38,color=c,edgecolors='white',linewidth=.5,zorder=3)
 ax.annotate(row['id'],(row['recovery_pct'],row['wrong']),xytext=off,textcoords='offset points',fontsize=8,fontweight='bold')
ax.set_xlabel('Exact 110-nt bodies / all designs (%)');ax.set_ylabel('Incorrect outputs');ax.set_xlim(95.05,98.4);ax.set_ylim(-10,190);ax.set_yticks([0,50,100,150]);clean(ax)
ax=fig.add_subplot(gs[1,1]);title(ax,'d','Complete-file validation');diagram(ax)
ax.text(.0,.93,'Input data set',fontweight='bold');ax.text(.57,.93,'File outcome',fontweight='bold')
for yy,name,status,color in zip([.74,.56,.38,.20,.02],['DNA Fountain','Photo File 1','Photo File 2','Photo File 3','Organick'],['MD5 matched','Bytes + ZIP CRCs','Bytes + ZIP CRCs','Incomplete','Not tested'],[GREEN,GREEN,GREEN,RED,DARK]):
 ax.text(.0,yy,name,va='center',fontsize=8);ax.text(.57,yy,status,va='center',fontsize=8,color=color,fontweight='medium')
plt.close(fig)

# S3: a sequence layout, index-level overlap, and the independently verified archive endpoint.
fig=plt.figure(figsize=(7.35,4.7));gs=fig.add_gridspec(2,2,height_ratios=[1,1.3],hspace=.6,wspace=.65)
ax=fig.add_subplot(gs[0,:]);title(ax,'a','DNA Fountain read layout and seeding');ax.set_xlim(-18,158);ax.set_ylim(-.3,3.7);ax.axis('off')
for start,length,label,color in [(0,16,'Seed 16 nt',GREEN),(16,128,'Fountain payload 128 nt',BLUE),(144,8,'RS 8 nt',PURPLE)]:
 ax.add_patch(Rectangle((start,2.35),length,.5,facecolor=color,edgecolor='white'))
 if length>20:ax.text(start+length/2,2.6,label,color='white',ha='center',va='center')
 else:ax.annotate(label,(start+length/2,2.9),xytext=(start+length/2,3.3),ha='center',fontsize=8,arrowprops={'arrowstyle':'-','lw':.5})
ax.add_patch(Rectangle((-10,1.25),10,.4,facecolor=GRAY));ax.text(-14,1.45,'10-nt\nanchor',ha='center',va='center',fontsize=7.4)
ax.annotate('',(151,1.45),(0,1.45),arrowprops={'arrowstyle':'->','color':ORANGE,'lw':2});ax.text(70,1.75,'R1  positions 0–150',ha='center')
ax.annotate('',(152,.55),(1,.55),arrowprops={'arrowstyle':'->','color':BLUE,'lw':2});ax.text(76,.82,'Reverse complement of R2  positions 1–151',ha='center')
ax.text(76,-.05,'k = 21, c = 0, r = 6; detect-only RS validation',ha='center',fontsize=8)
ax=fig.add_subplot(gs[1,0]);title(ax,'b','Overlap of output indices')
vals=[D['fountain'][k] for k in ['shared_indices','DBGPS2_only','merge_only']];ys=np.arange(3)
ax.barh(ys,vals,color=[GRAY,BLUE,ORANGE],height=.55);ax.set_yticks(ys,['Both methods','DBGPS2 only','Read merging only']);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlim(1,150000);ax.set_xlabel('Distinct indices (log scale)')
for y,v in zip(ys,vals):ax.text(v*1.14,y,f'{v:,}',va='center',fontsize=8)
clean(ax)
ax=fig.add_subplot(gs[1,1]);title(ax,'c','Archive validation');diagram(ax)
for y,t in [(.9,'71,990 exact graph outputs'),(.53,'67,088 fountain segments solved'),(.16,'2,146,816-byte archive')]:ax.text(.5,y,t,ha='center',va='center',fontsize=8.6,color=DARK)
for y1,y2 in [(.8,.64),(.43,.27)]:ax.annotate('',(.5,y2),(.5,y1),arrowprops={'arrowstyle':'->','color':GREEN,'lw':1.3})
ax.text(.5,-.06,'Published MD5 matched',ha='center',color=GREEN,fontweight='bold',fontsize=8.7)
save(fig,'FigS3_DNAFountain_validation')

# S4: candidate coverage and actual outer-code budgets, with explicit reference assistance.
fig=plt.figure(figsize=(7.35,6.5));gs=fig.add_gridspec(3,2,height_ratios=[.65,1.6,1.2],hspace=.8,wspace=.68)
ax=fig.add_subplot(gs[0,:]);title(ax,'a','Photolithographic code and reconstruction decisions');diagram(ax)
strand(ax,.01,.58,60,[(42,'Payload 42 nt',BLUE),(12,'Index 12 nt',GREEN),(6,'RS',PURPLE)],width=.86,height=.23)
ax.text(.91,.70,'60 nt',va='center',fontsize=8)
for xx,t in zip([.01,.35,.71],['Graph candidates','One output per index','Outer RS + archive']):ax.text(xx,.13,t,fontsize=8)
for x1,x2 in [(.24,.32),(.61,.68)]:ax.annotate('',(x2,.18),(x1,.18),arrowprops={'arrowstyle':'->','lw':.9})
ax=fig.add_subplot(gs[1,0]);title(ax,'b','Block 12 correct-candidate unions')
vals=[x['correct_candidates'] for x in D['block12_candidates']];ys=np.arange(5);labels=['Original pool\n179 settings','Assigned re-pool\nmatched 10','Separate pools\nmatched 10','Separate pool\n179 settings','Separate pool\n353 settings']
ax.barh(ys,vals,color=[GRAY,ORANGE,BLUE,BLUE,BLUE],height=.58)
ax.set_yticks(ys,labels,fontsize=7.2);ax.invert_yaxis();ax.set_xlim(0,12500);ax.set_xlabel('Indices with ≥1 exact candidate');ax.axvline(10977,color=RED,ls='--',lw=1)
for y,v in zip(ys,vals):ax.text(v-200,y,f'{v:,}',va='center',ha='right',color=DARK if y==0 else 'white',fontsize=7.5)
ax.text(10977,-.64,'10,977 needed',ha='right',color=RED,fontsize=7.2);clean(ax)
ax=fig.add_subplot(gs[1,1]);title(ax,'c','Actual selected block outputs')
bl=D['photo_selected_blocks'];xx=np.arange(1,13);yy=[r['budget'] for r in bl]
ax.bar(xx,yy,color=[BLUE]*11+[RED],width=.72);ax.axhline(5406,color=RED,ls='--',lw=1);ax.text(.7,5550,'RS limit 5,406',color=RED,fontsize=7.5)
ax.text(12,6600,'6,425',ha='center',fontsize=7.7,color=RED);ax.set_ylim(0,7400);ax.set_xlim(.4,12.6);ax.set_xticks([1,3,5,7,9,12]);ax.set_xlabel('Outer-code block');ax.set_ylabel(r'Worst-lane $2E_j + S$');clean(ax)
ax=fig.add_subplot(gs[2,:]);title(ax,'d','Block 12 supplementation controls')
controlvals=[b12['max_lane_budget'] for b12 in [D['block12_controls'][0]['after'],D['block12_controls'][1]['after']]]
values=[6425]+controlvals;labels=['Actual output\n9,973 exact; 15 wrong','Actual + 1,019\ndesign strands','Truth union + 318\ndesign strands']
for i,(v,label) in enumerate(zip(values,labels)):
 ax.barh(i,v,color=RED if i==0 else GRAY,edgecolor=RED if i==0 else DARK,hatch=None if i==0 else '///',height=.52,lw=.55)
 ax.text(v+70,i,f'{v:,}',va='center',fontsize=8)
ax.set_yticks(range(3),labels,fontsize=7.6);ax.invert_yaxis();ax.set_xlim(0,7300);ax.set_xlabel(r'Worst-lane $2E_j + S$');ax.axvline(5406,color=RED,ls='--',lw=1)
ax.text(5406,-.65,'5,406',color=RED,ha='center',fontsize=7.5);clean(ax)
save(fig,'FigS4_photolithographic_decoding')

# S5: actual operating points and truth-dependent graph diagnostics.
fig=plt.figure(figsize=(7.35,5.7));gs=fig.add_gridspec(3,2,height_ratios=[.8,1.4,1.35],hspace=.85,wspace=.65)
ax=fig.add_subplot(gs[0,:]);title(ax,'a','Organick strand format and validity check');diagram(ax)
strand(ax,.01,.56,150,[(20,'Primer',GRAY),(110,'Rotating ternary body 110 nt',BLUE),(20,'Primer',GRAY)],width=.97,height=.28)
ax.text(.5,.14,'Bidirectional graph search → format predicate → one distinct body',ha='center',fontsize=8.3)
ax=fig.add_subplot(gs[1:,0]);title(ax,'b','Reconstruction operating points')
rows=D['organick'];ys=np.arange(6)
for y,row,col in zip(ys,rows,colors):
 ax.scatter(row['recovery_pct'],y,color=col,s=36,zorder=3,edgecolors=DARK,lw=.3)
 ax.text(98.98,y,f"{row['wrong']:,}",va='center',ha='right',fontsize=7.7)
ax.text(98.98,-.9,'Wrong',ha='right',fontsize=7.6,fontweight='bold')
ax.set_yticks(ys,['A Initial\nforward','B Tuned\nbidirectional','C Inner code\nsingle run','D Ensemble\nmax recovery','E Ensemble\nhigh precision','F Ensemble\nzero errors'],fontsize=7.4);ax.invert_yaxis();ax.set_xlim(95,99.3);ax.set_ylim(5.6,-1.05);ax.set_xticks([95,96,97,98]);ax.set_xlabel('Exact bodies / all designs (%)');ax.grid(axis='x',color='#e5e7eb',lw=.5);clean(ax)
ax=fig.add_subplot(gs[1:,1]);title(ax,'c','Restricted path-existence ceiling')
bs=D['organick_path_bound']['bounds'];cx=np.array([r['cutoff']+1 for r in bs]);cy=np.array([r['either_k_either_direction']/607150*100 for r in bs])
ax.plot(cx,cy,'o-',color=DARK,ms=4,label='Truth-path diagnostic')
ax.scatter([1],[rows[3]['recovery_pct']],marker='*',color=BLUE,s=90,zorder=4,label='Best exact output')
ax.annotate('596,189 true paths\n596,006 exact outputs\nGap: 183 strands',(1,cy[0]),xytext=(1.7,96.5),fontsize=7.4,va='top',arrowprops={'arrowstyle':'-','lw':.65,'color':DARK})
ax.set_xlim(.7,7.3);ax.set_ylim(94,98.65);ax.set_xticks(range(1,8));ax.set_xlabel('Minimum k-mer count');ax.set_ylabel('Designed bodies with a true path (%)');ax.legend(frameon=False,loc='lower left',fontsize=7.2,handlelength=1.1);clean(ax)
save(fig,'FigS5_Organick_operating_points')
print('Generated',len(list(OUT.glob('*.png'))),'figures in PNG, PDF and SVG')
