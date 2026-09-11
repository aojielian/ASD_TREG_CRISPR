#!/usr/bin/env python3
"""Render the four final main figures from public source tables.

The renderer does not run inferential statistics or define new thresholds,
gene sets, effect estimates, or datasets. Counts, proportions, sorting and
empirical CDF coordinates are display operations.
"""
import argparse
import textwrap
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from figure_style import *

def truth(s):
    return s.astype(str).str.upper().isin(['TRUE', 'YES', '1'])

def numeric(s):
    return pd.to_numeric(s, errors='raise')

def source(df, name):
    match = df['source_file'].astype(str).map(lambda x: Path(x).name) == name
    out = df.loc[match].copy()
    if out.empty:
        raise ValueError(f'Missing frozen input in figure composite: {name}')
    return out

def records(df, record_type):
    out=df.loc[df.record_type==record_type].copy()
    if out.empty:raise ValueError(f'Missing record_type: {record_type}')
    return out

def orientation(ax, blocks):
    blank(ax)
    n=len(blocks)
    for i,(headline,detail,color) in enumerate(blocks):
        x=i/n
        ax.plot([x,x+.012],[.1,.82],color=color,lw=2.4)
        ax.text(x+.027,.75,headline,fontsize=8.5,fontweight='bold',va='top',color=color)
        ax.text(x+.027,.38,detail,fontsize=7.5,va='top',color=DARK,linespacing=1.4)
        if i<n-1:
            ax.annotate('',xy=((i+1)/n-.024,.48),xytext=((i+1)/n-.067,.48),
                        arrowprops={'arrowstyle':'->','color':GRAY,'lw':.8})

def enrichment(ax, old, color, xmax):
    props=old.loc[old.row_type=='proportion']
    effect=old.loc[old.row_type=='effect'].iloc[0]
    blank(ax)
    bar_left,bar_width=.02,.72
    for j,(_,r) in enumerate(props.iterrows()):
        y=.81-j*.27
        label='ASD genes' if 'SFARI' in r['label'] else 'Other genes'
        num,den=int(float(r.numerator)),int(float(r.denominator))
        value=float(r.y)*100
        ax.text(.02,y+.10,f'{label}  {num:,} / {den:,}',fontsize=8,
                weight='bold' if j==0 else 'normal')
        ax.add_patch(Rectangle((bar_left,y-.035),bar_width*value/xmax,.065,
                               color=color if j==0 else GRAY,lw=0))
        ax.text(.80,y,f'{value:.1f}%',fontsize=8,va='center')
    ax.plot([bar_left,bar_left+bar_width],[.405,.405],color=DARK,lw=.55)
    ticks=[0,3,6] if xmax==7 else [0,5,10,15]
    for tick in ticks:
        xx=bar_left+bar_width*tick/xmax
        ax.plot([xx,xx],[.405,.38],color=DARK,lw=.55)
        ax.text(xx,.35,str(tick),ha='center',va='top',fontsize=7)
    ax.text(.80,.405,'Hits (%)',fontsize=7,va='center')
    ax.text(.02,.15,f'OR {float(effect.effect):.2f}',fontsize=10,weight='bold',color=color)
    ax.text(.02,-.01,'One-sided P = '+pformat(effect.p,4),fontsize=7.5)

def figure1(df,out):
    full=records(df,'screen_gene')
    old=records(df,'enrichment')
    full['author_lfc']=numeric(full.author_lfc);full['screen_strength']=numeric(full.screen_strength)
    asd=truth(full.sfari_score1); hit=truth(full.screen_hit)
    regs=full.loc[asd & hit].copy().sort_values('author_lfc')
    assert len(full)==18966 and asd.sum()==242 and len(regs)==14
    fig=new_figure('ASD risk genes converge on FOXP3 regulation',6.6)
    a=panel(fig,[.075,.825,.875,.077],'a','Published screen integration',1.10)
    orientation(a,[('ASD risk genes','SFARI score 1\n242 measured genes',NAVY),
                   ('Human FOXP3 screen','Genome-wide CRISPR\n18,966 measured genes',DARK),
                   ('Functional regulators','14 ASD genes\nBoth effect directions',TEAL)])
    b=panel(fig,[.09,.135,.47,.60],'b','ASD genes within the full screen')
    clean(b)
    b.scatter(full.loc[~asd,'author_lfc'],full.loc[~asd,'screen_strength'],s=4,
              c=GRAY,alpha=.42,linewidths=0)
    b.scatter(full.loc[asd & ~hit,'author_lfc'],full.loc[asd & ~hit,'screen_strength'],
              s=10,facecolor='white',edgecolor=NAVY,linewidth=.45,alpha=.7)
    for direction,color in [('SUPPORTS_FOXP3',TEAL),('RESTRAINS_FOXP3',ORANGE)]:
        r=regs.loc[regs.regulator_direction==direction]
        b.scatter(r.author_lfc,r.screen_strength,s=23,color=color,edgecolor='white',lw=.4,zorder=5)
    b.axvline(0,color=GRAY,lw=.7,ls=(0,(3,3)))
    b.set_xlabel('Author screen effect (LFC)')
    b.set_ylabel(r'Screen strength: $-\log_{10}(P_{\mathrm{min\ tail}})$')
    b.set_ylim(-.12,float(full.screen_strength.max())+.9)
    b.set_xlim(float(full.author_lfc.min())-.2,float(full.author_lfc.max())+.3)
    # Labels are fixed concept-selected genes; offsets do not change plotted coordinates.
    offsets={'SMARCA4':(9,11),'EP300':(16,5),'FOXP1':(13,-11),
             'TBL1XR1':(0,19),'CREBBP':(-14,-22)}
    for gene,(dx,dy) in offsets.items():
        r=regs.loc[regs.gene==gene].iloc[0]
        color=TEAL if r.regulator_direction=='SUPPORTS_FOXP3' else ORANGE
        b.annotate(gene,(r.author_lfc,r.screen_strength),xytext=(dx,dy),textcoords='offset points',
                   fontsize=7,fontstyle='italic',color=color,weight='bold',
                   arrowprops={'arrowstyle':'-','color':color,'lw':.5},zorder=6)
    leg=marker_legend(b,[('Other genes',GRAY,'o'),('ASD, not a hit',NAVY,'o'),
                     ('ASD, supports FOXP3',TEAL,'o'),('ASD, restrains FOXP3',ORANGE,'o')],
                  loc='upper left',bbox_to_anchor=(-.01,-.13),ncol=2,columnspacing=.8)
    leg.legend_handles[1].set_markerfacecolor('white')
    leg.legend_handles[1].set_markeredgewidth(.6)
    c=panel(fig,[.685,.57,.28,.17],'c','Regulator enrichment',1.14)
    enrichment(c,old,NAVY,7)
    d=panel(fig,[.755,.135,.20,.34],'d','ASD regulator effects',1.14,header_left=.685)
    clean(d,grid='x')
    for i,(_,r) in enumerate(regs.iterrows()):
        col=TEAL if r.regulator_direction=='SUPPORTS_FOXP3' else ORANGE
        d.hlines(i,0,r.author_lfc,color=col,lw=1)
        d.scatter(r.author_lfc,i,color=col,s=17,zorder=4)
    d.axvline(0,color=GRAY,lw=.6)
    d.set_yticks(range(len(regs)),regs.gene,fontsize=7,fontstyle='italic')
    d.tick_params(axis='y',length=0,pad=4)
    d.invert_yaxis();d.set_ylim(len(regs)-.5,-.5)
    d.set_xlim(-1.55,2.05);d.set_xticks([-1,0,1,2]);d.set_xlabel('Author screen effect (LFC)',fontsize=7)
    d.text(0,-.27,'Restrains',transform=d.transAxes,ha='left',color=ORANGE,fontsize=7)
    d.text(1,-.27,'Supports',transform=d.transAxes,ha='right',color=TEAL,fontsize=7)
    save(fig,out,1)

def figure2(df,out):
    old=records(df,'enrichment')
    members=records(df,'common_asd_gene')
    assert truth(members.chen_measured).all() and truth(members.umhoefer_treg_measured).all()
    members['label']=members.gene
    members['group']=members.cross_screen_category
    members['x']=members.chen_hit
    members['y']=members.umhoefer_treg_hit
    order={'BOTH_HIT':0,'UMHOEFER_TREG_ONLY':1,'NEITHER':2}
    members['_order']=members.group.map(order)
    members=members.sort_values(['_order','label']).reset_index(drop=True)
    assert len(members)==33 and (members.group=='BOTH_HIT').sum()==2
    fig=new_figure('Independent support in primary human Treg cells',5.9)
    a=panel(fig,[.075,.825,.875,.075],'a','Independent screen context',1.10)
    orientation(a,[('Genome-wide discovery','Human induced Tregs\nFOXP3 phenotype',NAVY),
                   ('Common ASD universe','33 genes tested\nin both screens',DARK),
                   ('Primary-Treg screen','Targeted CRISPR\n1,349 measured genes',TEAL)])
    b=panel(fig,[.15,.505,.805,.225],'b','Hit membership across 33 shared ASD genes',1.11,header_left=.075)
    b.set_xlim(-.6,32.6);b.set_ylim(-.68,2.55);b.axis('off')
    for start,end,col,label in [(0,1,NAVY,'Both: 2'),(2,4,TEAL,'Primary only: 3'),(5,32,GRAY,'Neither: 28')]:
        b.plot([start-.35,end+.35],[2.33,2.33],color=col,lw=2)
        # Separate headers for first two small groups avoids overlapping group names.
        xpos=(start+end)/2
        # Group boundary bars follow the exact ordered category memberships.
    for i,r in members.iterrows():
        b.text(i,-.24,r.label,rotation=90,fontsize=7,fontstyle='italic',ha='center',va='top',
               color=NAVY if i<5 else DARK)
        for y,col,val in [(1.5,NAVY,r.x),(.6,TEAL,r.y)]:
            b.scatter(i,y,s=37,marker='s',facecolor=col if val=='YES' else 'white',
                      edgecolor=col if val=='YES' else GRAY,linewidth=.65)
    for y,label,col in [(1.5,'Discovery',NAVY),(.6,'Primary Treg',TEAL)]:
        b.text(-1.05,y,label,ha='right',va='center',fontsize=7.5,color=col)
    b.text(0,2.65,'2 both-hit  ·  3 primary-Treg-only  ·  28 neither',fontsize=7.5,va='center')
    c=panel(fig,[.085,.115,.35,.235],'c','Independent enrichment',1.11)
    enrichment(c,old,TEAL,18)
    d=panel(fig,[.57,.11,.385,.24],'d','Key-gene membership',1.11)
    blank(d)
    xs=[.37,.54,.75,.92]
    d.text(.455,.98,'Discovery',ha='center',weight='bold',fontsize=7.5,color=NAVY)
    d.text(.835,.98,'Primary Treg',ha='center',weight='bold',fontsize=7.5,color=TEAL)
    for x,label in zip(xs,['Tested','Hit','Tested','Hit']):d.text(x,.86,label,ha='center',fontsize=7)
    lookup=members.set_index('label')
    for i,gene in enumerate(['FOXP1','ZMYND8','CTCF','PTEN','SATB1']):
        y=.73-i*.14
        d.text(.0,y,gene,va='center',fontsize=8,fontstyle='italic')
        for j,(x,value) in enumerate(zip(xs,[True,lookup.loc[gene,'x']=='YES',True,lookup.loc[gene,'y']=='YES'])):
            color=NAVY if j<2 else TEAL
            d.scatter(x,y,s=25,marker='o',facecolor=color if value else 'white',edgecolor=color,lw=.7)
        table_rule(d,y-.069)
    d.text(0,-.10,'Filled = yes     Open = no',fontsize=7,color=DARK)
    save(fig,out,2)

def membership_map(ax,common):
    blank(ax)
    ax.set_xlim(-.19,2.02);ax.set_ylim(-.30,2.11)
    ch=truth(common.CHEN_HIT_BOOL);sh=truth(common.SHIFRUT_HIT);asd=truth(common.SFARI_SCORE1)
    for c in [0,1]:
        for s in [0,1]:
            sub=common.loc[(ch==bool(c)) & (sh==bool(s))]
            asub=sub.loc[truth(sub.SFARI_SCORE1)]
            color=TEAL if (c,s)==(1,0) else ORANGE if (c,s)==(0,1) else NAVY
            ax.add_patch(Rectangle((c,s),.97,.97,facecolor=PALE,edgecolor='white',lw=1))
            ax.plot([c+.025,c+.945],[s+.925,s+.925],color=color,lw=1.3)
            label={(1,0):'FOXP3 only',(0,1):'TCR only',(1,1):'Both hits',(0,0):'Neither hit'}[(c,s)]
            ax.text(c+.06,s+.84,label,fontsize=8,weight='bold',color=color)
            ax.text(c+.06,s+.72,f'{len(asub)} ASD / {len(sub):,} all',fontsize=7,color=DARK)
            if len(asub)<=12:
                genes=sorted(asub.approved_symbol)
                for i,g in enumerate(genes):
                    col=i//6; row=i%6
                    ax.text(c+.06+col*.46,s+.57-row*.09,g,fontsize=7,fontstyle='italic',color=color)
            else:
                n=len(asub);cols=20
                xx=c+.065+(np.arange(n)%cols)*.042
                yy=s+.56-(np.arange(n)//cols)*.039
                ax.scatter(xx,yy,s=2.9,c=GRAY,marker='s',linewidths=0)
                ax.text(c+.06,s+.045,'1 square = 1 ASD gene',fontsize=7,color=DARK)
    ax.text(.485,-.11,'No hit',ha='center',fontsize=7)
    ax.text(1.485,-.11,'Hit',ha='center',fontsize=7)
    ax.text(.985,-.25,'FOXP3 screen (Chen)',ha='center',fontsize=8)
    ax.text(-.05,.485,'No hit',rotation=90,va='center',ha='right',fontsize=7)
    ax.text(-.05,1.485,'Hit',rotation=90,va='center',ha='right',fontsize=7)
    ax.text(-.18,.985,'TCR / proliferation (Shifrut)',rotation=90,va='center',ha='right',fontsize=8)

def figure3(df,out):
    common=records(df,'common_gene')
    direct=records(df,'preference_summary').iloc[0]
    classes=records(df,'class_conditioned').set_index('screen')
    joint=records(df,'joint_conditioned').set_index('screen')
    assert len(common)==18615 and truth(common.SFARI_SCORE1).sum()==242
    fig=new_figure('FOXP3 phenotype preference and robustness',6.45)
    a=panel(fig,[.075,.825,.875,.076],'a','Direct comparison within a common universe',1.10)
    orientation(a,[('FOXP3 regulation','Chen human screen',TEAL),
                   ('Common tested universe','18,615 genes · 242 ASD',DARK),
                   ('TCR / proliferation','Shifrut human screen',ORANGE)])
    b=panel(fig,[.09,.145,.45,.585],'b','Observed cross-phenotype membership',1.03)
    membership_map(b,common)
    c=panel(fig,[.64,.485,.325,.245],'c','Phenotype-preference categories',1.08)
    blank(c)
    c.text(.0,.91,'Category',fontsize=7,color=DARK)
    c.text(.63,.91,'ASD',fontsize=7,color=NAVY,ha='right',weight='bold')
    c.text(.98,.91,'All genes',fontsize=7,color=DARK,ha='right')
    cats=[('FOXP3_ONLY','FOXP3-only',TEAL,'d_plus'),('CONCORDANT','Concordant',NAVY,'d_zero'),('TCR_ONLY','TCR-only',ORANGE,'d_minus')]
    for i,(key,label,col,field) in enumerate(cats):
        y=.74-i*.18
        # Category spelling is validated by numeric D; no threshold is introduced.
        dv={'d_plus':'1','d_zero':'0','d_minus':'-1'}[field]
        n_asd=int((truth(common.SFARI_SCORE1) & (common.D==dv)).sum())
        n_all=int(float(direct[field]))
        c.text(.0,y,label,fontsize=8,color=col)
        c.text(.63,y,f'{n_asd:,}',ha='right',fontsize=10,color=col,weight='bold')
        c.text(.98,y,f'{n_all:,}',ha='right',fontsize=8,color=DARK)
        table_rule(c,y-.055)
    c.text(0,.08,'Exact preference P = '+pformat(direct.exact_null_p,4),fontsize=7.5)
    c.text(0,-.03,'Concordant = both-hit + neither-hit',fontsize=6.9,color=DARK)
    d=panel(fig,[.64,.105,.325,.27],'d','Conditioning and sensitivity',1.10)
    blank(d)
    d.text(0,.93,'Regulator-class conditioned',fontsize=8,weight='bold',color=TEAL)
    d.text(.56,.79,'Exact P',fontsize=7,ha='center');d.text(.92,.79,'Holm P',fontsize=7,ha='center')
    for i,(screen,label) in enumerate([('CHEN','Chen'),('UMHOEFER_TREG','Primary Treg')]):
        r=classes.loc[screen];y=.66-i*.12
        d.text(0,y,label,fontsize=7.5)
        d.text(.56,y,f'{float(r.EXACT_ONE_SIDED_P):.8f}',ha='center',fontsize=7)
        d.text(.92,y,f'{float(r.holm_p_two_screen):.7f}',ha='center',fontsize=7)
    table_rule(d,.47)
    d.text(0,.37,'Strict joint-conditioned',fontsize=8,weight='bold',color=DARK)
    d.text(.56,.25,'Exact P',fontsize=7,ha='center');d.text(.92,.25,'Coverage',fontsize=7,ha='center')
    for i,(screen,label) in enumerate([('CHEN','Chen'),('UMHOEFER_TREG','Primary Treg')]):
        r=joint.loc[screen];y=.12-i*.12
        d.text(0,y,label,fontsize=7.5)
        d.text(.56,y,f'{float(r.EXACT_ONE_SIDED_P):.7f}',ha='center',fontsize=7)
        d.text(.92,y,f'{float(r.COMPLETE_CASE_ASD_COVERAGE)*100:.1f}%',ha='center',fontsize=7)
    d.text(0,-.14,'Restrictive sensitivity / coverage boundary',fontsize=6.8,color=DARK)
    save(fig,out,3)

def figure4(df,out):
    summary=records(df,'program')
    genes=records(df,'patient_gene')
    summary['_id']=summary.program_id.str.extract(r'(\d+)').astype(int)
    summary=summary.sort_values('_id').set_index('_id')
    effects=numeric(genes.signed_effect)
    assert len(genes)==5851
    supported=summary.index[summary.interpretation_status=='FORMAL_SUPPORTED_Q_LT_005'].tolist()
    assert supported==[1,6,8,9]
    fig=new_figure('Patient Treg transcriptomic context',7.05)
    a=panel(fig,[.065,.35,.53,.55],'a','Program evidence atlas',1.025)
    blank(a)
    a.set_xlim(0,1);a.set_ylim(-.50,10.8)
    a.text(0,10.3,'Program / source descriptor',fontsize=7,color=DARK)
    a.text(.62,10.3,'Coverage',fontsize=7,color=DARK,ha='center')
    a.text(.83,10.3,'Direction',fontsize=7,color=DARK,ha='center')
    a.text(.96,10.3,'Status',fontsize=7,color=DARK,ha='center')
    for i,r in summary.iterrows():
        y=10-i
        sup=i in supported;eligible=r.formal_eligibility=='ELIGIBLE'
        col=ORANGE if sup and r.effect_direction=='ASD_HIGHER' else TEAL if sup else NAVY if eligible else GRAY
        label=textwrap.fill(r.program_short_label,width=22,break_long_words=False)
        a.text(0,y+.14,f'M{i}',fontsize=8.5,weight='bold',color=NAVY,va='center')
        a.text(.105,y+.14,label,fontsize=7,color=DARK,va='center',linespacing=1.25)
        a.add_patch(Rectangle((.52,y+.07),.20,.17,facecolor=LIGHT,lw=0))
        a.add_patch(Rectangle((.52,y+.07),.20*float(r.coverage),.17,facecolor=col,lw=0,
                               hatch='////' if not eligible else None,edgecolor='white'))
        a.text(.62,y-.16,f'{int(float(r.in_valid_universe))}/{int(float(r.source_size))}',fontsize=6.9,ha='center',va='center')
        direction='↑' if r.effect_direction=='ASD_HIGHER' else '↓'
        a.text(.83,y+.14,direction,fontsize=12,ha='center',va='center',color=col)
        a.scatter(.96,y+.14,s=31,marker='o' if eligible else 's',
                  facecolor=col if sup or not eligible else 'white',edgecolor=col,lw=.8)
        table_rule(a,y-.43)
    # Status markers occupy their own, uncrowded line below the atlas.
    a.scatter([.015,.285],[-.70]*2,s=24,marker='o',c=[TEAL,'white'],edgecolors=[TEAL,NAVY],clip_on=False)
    a.scatter([.71],[-.70],s=24,marker='s',c=[GRAY],edgecolors=[GRAY],clip_on=False)
    for x,text_ in [(.045,'Supported'),(.315,'Eligible, unsupported'),(.74,'Ineligible')]:
        a.text(x,-.70,text_,fontsize=7,va='center')
    b=panel(fig,[.675,.35,.28,.55],'b','Observed gene distributions',1.025)
    blank(b)
    b.text(0,1.0,'Empirical cumulative fraction',fontsize=7,color=DARK,va='top')
    for j,i in enumerate(supported):
        y=.35+(3-j)*.135
        ax=fig.add_axes([.70,y,.255,.105])
        r=summary.loc[i];m=truth(genes[f'module_{i}_member'])
        col=ORANGE if r.effect_direction=='ASD_HIGHER' else TEAL
        assert int(m.sum())==int(float(r.in_valid_universe))
        ecdf(ax,effects[~m],GRAY,'Other measured genes',1)
        ecdf(ax,effects[m],col,'Program genes',1.15)
        ax.set_xscale('symlog',linthresh=2,linscale=1)
        ax.set_xlim(-4.2,15);ax.set_ylim(0,1.03)
        ax.set_xticks([-4,-1,1,4,10]);ax.set_xticklabels(['−4','−1','1','4','10'])
        ax.set_yticks([0,.5,1]);ax.set_yticklabels(['0','.5','1'])
        clean(ax,grid='y');ax.tick_params(labelsize=7,pad=1.5)
        ax.axvline(0,color=LIGHT,lw=.5,zorder=0)
        ax.text(.04,.90,f'M{i}',transform=ax.transAxes,color=col,fontsize=8,weight='bold',va='top')
        ax.text(.98,.12,'q = '+pformat(r.formal_q_BH10,3),transform=ax.transAxes,ha='right',fontsize=7,color=DARK)
        if j<3:ax.set_xticklabels([])
        else:ax.set_xlabel('Author signed fold change\n(symmetric-log axis)',fontsize=7,labelpad=2)
    for x0,x1,color in [(0,.04,ORANGE),(.04,.08,TEAL),(.55,.63,GRAY)]:
        b.plot([x0,x1],[-.15,-.15],color=color,lw=1.4,clip_on=False)
    b.text(.11,-.15,'Program',fontsize=7,color=DARK,va='center')
    b.text(.66,-.15,'Other features',fontsize=7,color=DARK,va='center')
    c=panel(fig,[.065,.075,.89,.125],'c','Perturbation programs in patient context',1.12)
    blank(c)
    for x,head,detail,col in [(.0,'FOXP3 perturbation','10 published programs',NAVY),
                            (.34,'Patient Treg transcriptome','5 coverage-eligible programs',DARK),
                            (.70,'4 supported programs','M1 ↑   M6 ↓   M8 ↓   M9 ↓',TEAL)]:
        c.text(x,.78,head,fontsize=8,weight='bold',color=col)
        c.text(x,.48,detail,fontsize=7.5,color=DARK)
    c.plot([.315,.315],[.42,.85],color=GRAY,lw=.6)
    c.plot([.675,.675],[.42,.85],color=GRAY,lw=.6)
    c.text(0,.10,'Program-level patient context complements the functional screens.',fontsize=8,color=NAVY)
    save(fig,out,4)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',required=True,type=Path)
    ap.add_argument('--output-dir',required=True,type=Path)
    ap.add_argument('--figures',default='1,2,3,4')
    args=ap.parse_args();configure()
    for i in map(int,args.figures.split(',')):
        df=pd.read_csv(args.input_dir/f'Figure{i}_source_data.tsv',sep='\t',dtype=str,keep_default_na=False)
        globals()[f'figure{i}'](df,args.output_dir)
        print(f'RENDERED=Figure{i}_FINAL.pdf,png,svg')

if __name__=='__main__':main()
