#!/usr/bin/env python3
"""Render Supplementary Figures S1-S5 from public source tables."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from supp_style import *

def records(df,name):
    out=df.loc[df.record_type==name].copy()
    if out.empty:raise ValueError(f'Missing record_type: {name}')
    return out
def truth(s):return s.astype(str).str.upper().isin(['TRUE','YES','1'])
def num(s):return pd.to_numeric(s,errors='raise')
def intval(x):return int(float(x))
def modules(df):
    df=df.copy();df['_m']=df.program_id.str.extract(r'(\d+)').astype(int)
    return df.sort_values('_m')
def title(n,text,height):return new_figure(f'Supplementary Figure S{n} | {text}',height)
def subtitle(ax,text,color=NAVY):ax.set_title(text,fontsize=8.5,color=color,pad=8)

def fraction_bars(ax,groups,xmax,xlabel='Hits (%)',color=TEAL):
    """Groups: label,numerator,denominator; fractions are display arithmetic only."""
    blank(ax)
    for j,(label,k,n) in enumerate(groups):
        y=.80-j*.28;value=100*k/n
        ax.text(0,y+.105,f'{label}   {k:,} / {n:,}',fontsize=8,
                weight='bold' if j==0 else 'normal')
        ax.add_patch(Rectangle((0,y-.025),.70*value/xmax,.065,
                               facecolor=color if j==0 else GRAY,lw=0))
        if k==0:ax.plot([0,0],[y-.025,y+.04],color=color,lw=1)
        ax.text(.76,y+.006,f'{value:.1f}%',fontsize=8,va='center')
    ax.plot([0,.7],[.37,.37],color=DARK,lw=.6)
    ticks=np.linspace(0,xmax,4)
    for tick in ticks:
        x=.7*tick/xmax;ax.plot([x,x],[.37,.34],color=DARK,lw=.6)
        ax.text(x,.31,f'{tick:g}',fontsize=7,ha='center',va='top')
    ax.text(.76,.37,xlabel,fontsize=7,va='center')

def figure_s1(df,out):
    genes=records(df,'common_asd_gene').sort_values('gene').reset_index(drop=True)
    er=records(df,'tconv_enrichment').iloc[0]
    assert len(genes)==33
    prefixes=['chen','umhoefer_treg','umhoefer_tconv']
    values=np.array([[float(r[p+'_author_lfc']) for p in prefixes] for _,r in genes.iterrows()])
    limit=float(np.abs(values).max())
    cmap=LinearSegmentedColormap.from_list('signed_lfc',[ORANGE,'#FFFFFF',TEAL])
    norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
    fig=title(1,'Cross-screen gene evidence',7.15)
    a=panel(fig,[.14,.13,.455,.715],'a','33 shared measured ASD genes',1.055,header_left=.075)
    a.set_xlim(-.05,6);a.set_ylim(32.6,-2.1);a.axis('off')
    for j,label in enumerate(['Chen','Primary Treg','Tconv']):
        a.text(j*2+.72,-1.78,label,ha='center',fontsize=8,weight='bold')
        a.text(j*2+.47,-.92,'LFC',ha='center',fontsize=7,color=DARK)
        a.text(j*2+1.44,-.92,'Hit',ha='center',fontsize=7,color=DARK)
    for i,r in genes.iterrows():
        a.text(-.19,i,r.gene,ha='right',va='center',fontsize=7,fontstyle='italic')
        for j,prefix in enumerate(prefixes):
            assert str(r[prefix+'_measured']).upper()=='TRUE'
            value=float(r[prefix+'_author_lfc'])
            cellcolor=cmap(norm(value))
            a.add_patch(Rectangle((j*2-.02,i-.43),1.02,.86,facecolor=cellcolor,edgecolor='white',lw=.3))
            display=f'{value:.2f}'
            if display=='-0.00':display='0.00'
            # Text contrast depends only on displayed RGB luminance, not science.
            rgb=np.array(cellcolor[:3]);linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
            luminance=float(linear @ np.array([.2126,.7152,.0722]))
            textcolor='white' if luminance<.28 else DARK
            a.text(j*2+.47,i,display,fontsize=6.9,ha='center',va='center',color=textcolor)
            hit=str(r[prefix+'_hit']).upper()=='YES'
            a.scatter(j*2+1.44,i,s=14,facecolor=NAVY if hit else 'white',edgecolor=NAVY,lw=.6)
    cbax=fig.add_axes([.22,.082,.27,.012])
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cbax,orientation='horizontal')
    cb.solids.set_rasterized(False)
    cb.set_ticks([-limit,0,limit]);cb.set_ticklabels([f'{-limit:.2f}','0',f'{limit:.2f}'])
    cb.ax.tick_params(labelsize=7,length=2,pad=2);cb.set_label('Signed author LFC',fontsize=7,labelpad=2)
    b=panel(fig,[.675,.555,.28,.28],'b','Assay context',1.14)
    blank(b)
    for j,(head,line,color) in enumerate([
        ('Chen','Genome-wide discovery\nHuman induced Tregs',NAVY),
        ('Primary Treg','Targeted screen\nIndependent Treg support',TEAL),
        ('Tconv','Targeted screen\nConventional T-cell context',DARK)]):
        y=.90-j*.32
        b.plot([0,0],[y-.18,y+.015],color=color,lw=2)
        b.text(.04,y,head,fontsize=8.5,weight='bold',color=color)
        b.text(.04,y-.055,line,fontsize=7.5,va='top',linespacing=1.3,color=DARK)
    b.scatter([.01,.01],[-.15,-.27],s=17,facecolors=[NAVY,'white'],edgecolors=NAVY,lw=.6,clip_on=False)
    b.text(.07,-.15,'Author hit',fontsize=7.5,va='center')
    b.text(.07,-.27,'Measured, not a hit',fontsize=7.5,va='center')
    c=panel(fig,[.675,.16,.28,.24],'c','Tconv enrichment',1.11)
    fraction_bars(c,[('ASD',intval(er.ASD_hits_or_observed),intval(er.ASD_n)),
                     ('Other',intval(er.comparator_hits_or_expected),intval(er.comparator_n))],15,color=DARK)
    c.text(0,.12,'One-sided Fisher P = '+pformat(er.P,3),fontsize=7.5)
    save(fig,out,1)

def figure_s2(df,out):
    recal=records(df,'recalibration');omission=records(df,'omission')
    fig=title(2,'Recalibration and influence sensitivity',6.8)
    host=panel(fig,[.075,.63,.875,.23],'a','Recalibration and removal of shared hits',1.15);blank(host)
    for j,(screen,label,color) in enumerate([('CHEN','Chen',NAVY),('UMHOEFER_TREG','Primary Treg',TEAL)]):
        x0=.11+j*.46
        ax=fig.add_axes([x0,.66,.36,.165]);clean(ax)
        r=recal.loc[recal.screen==screen].set_index('analysis_type_or_removal')
        selected=[r.loc['FULL_UNADJUSTED_RECALIBRATED'],r.loc['REMOVE_FOXP1;ZMYND8']]
        effects=[float(v.effect) for v in selected]
        ax.plot([0,1],effects,color=GRAY,lw=.7,zorder=1)
        ax.scatter([0,1],effects,s=25,c=color,zorder=3)
        ax.axhline(1,color=GRAY,lw=.7,ls=(0,(3,3)))
        ax.set_xlim(-.4,1.4);ax.set_ylim(0,10);ax.set_yticks([0,5,10]);ax.set_ylabel('Odds ratio')
        ax.set_xticks([0,1],['Recalibrated','Remove\nFOXP1 + ZMYND8'],fontsize=7)
        subtitle(ax,label,color)
        for k,v in enumerate(selected):
            ax.text(k,effects[k]+.6,f'{effects[k]:.2f}',ha='center',fontsize=7.5,color=color)
            ax.text(k,-.47,f'{intval(v.K)}/{intval(v.A)}\nP = {pformat(v.one_sided_P,3)}',
                    transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7,linespacing=1.4)
    host=panel(fig,[.075,.355,.875,.10],'b','Recalibrated class-conditioned enrichment',1.43);blank(host)
    for j,(screen,label,color) in enumerate([('CHEN','Chen',NAVY),('UMHOEFER_TREG','Primary Treg',TEAL)]):
        ax=fig.add_axes([.11+j*.46,.375,.36,.075]);clean(ax,grid='x')
        r=recal.loc[(recal.screen==screen)&(recal.analysis_type_or_removal=='CLASS_CONDITIONED_RECALIBRATED')].iloc[0]
        value=float(r.observed_expected_ratio)
        ax.axvline(1,color=GRAY,lw=.7,ls=(0,(3,3)))
        ax.scatter(value,.5,c=color,s=25,zorder=3)
        ax.text(value,.85,f'{value:.2f}',fontsize=8,ha='center',color=color)
        ax.text(0,1.17,label,transform=ax.transAxes,fontsize=8,color=color,weight='bold')
        ax.set_ylim(0,1.35);ax.set_yticks([]);ax.set_xlim(0,4);ax.set_xticks([0,1,2,3,4]);ax.set_xlabel('Observed / expected',fontsize=7)
    host=panel(fig,[.075,.09,.875,.145],'c','Single-gene omission sensitivity',1.38);blank(host)
    for j,(screen,label,color,expected) in enumerate([('CHEN','Chen',NAVY,242),('UMHOEFER_TREG','Primary Treg',TEAL,33)]):
        ax=fig.add_axes([.11+j*.46,.09,.36,.145]);clean(ax)
        sub=omission.loc[omission.screen==screen].sort_values('removed_gene').reset_index(drop=True)
        assert len(sub)==expected
        ishit=truth(sub.removed_was_author_hit)
        y=-np.log10(num(sub.one_sided_fisher_p))
        ax.scatter(np.arange(len(sub)),y,s=9,c=np.where(ishit,color,GRAY),lw=0)
        ax.axhline(-np.log10(.05),color=GRAY,lw=.7,ls=(0,(3,3)))
        ax.text(.02,.245,'P = 0.05',transform=ax.transAxes,fontsize=6.8,color=DARK)
        ax.set_ylim(0,6);ax.set_yticks([0,2,4,6]);ax.set_xlim(-2,len(sub)+1)
        ax.set_xticks([0,len(sub)-1],['1',str(len(sub))])
        ax.set_xlabel('Removed gene (alphabetical order)',fontsize=7)
        ax.set_ylabel(r'$-\log_{10}$(Fisher P)',fontsize=7)
        subtitle(ax,f'{label} · {expected} omissions',color)
    handles=[Line2D([],[],marker='o',ls='none',color=NAVY,markersize=4,label='Chen: removed hit'),
             Line2D([],[],marker='o',ls='none',color=TEAL,markersize=4,label='Primary Treg: removed hit'),
             Line2D([],[],marker='o',ls='none',color=GRAY,markersize=4,label='Removed nonhit')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.002),ncol=3,frameon=False,fontsize=7)
    save(fig,out,2)

def figure_s3(df,out):
    strata=records(df,'class_stratum');coverage=records(df,'joint_coverage');mouse=records(df,'mouse_sensitivity').set_index('row_type')
    fig=title(3,'Regulator classes and analysis scope',6.5)
    host=panel(fig,[.075,.50,.875,.34],'a','Regulator-class strata',1.16);blank(host)
    groups=[('CHEN','TRUE','Chen / library in',NAVY),('CHEN','FALSE','Chen / library out',NAVY),('UMHOEFER_TREG','TRUE','Primary Treg',TEAL)]
    classes=['BOTH_LISTED','EPIFACTOR_ONLY','TF_ONLY','NEITHER_LISTED']
    labels=['Both (TF + EpiFactors)','EpiFactors only','TF only','Neither']
    for j,(screen,library,label,color) in enumerate(groups):
        ax=fig.add_axes([.245+j*.247,.51,.205,.31]);clean(ax,grid='x')
        sub=strata.loc[(strata.screen==screen)&(strata.library_member==library)].set_index('class_label')
        for k,cl in enumerate(classes):
            r=sub.loc[cl];n,h,a,hit=[intval(r[v]) for v in ['N','H','A','K']]
            for offset,kk,nn,fill in [(-.16,hit,a,True),(.16,h-hit,n-a,False)]:
                value=100*kk/nn
                ax.barh(k+offset,value,height=.23,color=color if fill else 'white',edgecolor=color,lw=.7)
                if value==0:ax.scatter(0,k+offset,s=10,color=color if fill else 'white',edgecolor=color,lw=.6,zorder=4)
                ax.text(103,k+offset,f'{kk}/{nn:,}',fontsize=6.8,ha='left',va='center',color=DARK)
        ax.set_xlim(-2,145);ax.set_xticks([0,50,100]);ax.set_ylim(3.6,-.6)
        ax.set_yticks(range(4),labels if j==0 else ['']*4,fontsize=7)
        ax.tick_params(axis='y',length=0,pad=6);ax.set_xlabel('Hit rate (%)',fontsize=7)
        subtitle(ax,label,color)
    fig.text(.245,.443,'Filled: ASD genes     Open: other genes     Labels: hits / measured',fontsize=7,color=DARK)
    b=panel(fig,[.085,.12,.35,.225],'b','Joint-model complete cases',1.16);blank(b)
    for j,(screen,label,color) in enumerate([('CHEN','Chen',NAVY),('UMHOEFER_TREG','Primary Treg',TEAL)]):
        r=coverage.loc[coverage.screen==screen].iloc[0]
        k,n=intval(r.COMPLETE_CASE_ASD_MEASURED),intval(r.ORIGINAL_MEASURED_ASD)
        value=100*float(r.COMPLETE_CASE_ASD_COVERAGE);y=.81-j*.36
        b.text(0,y+.10,f'{label}   {k}/{n}',fontsize=8,color=color,weight='bold')
        b.add_patch(Rectangle((0,y-.015),.70*value/100,.07,facecolor=color,lw=0))
        b.text(.76,y+.02,f'{value:.1f}%',fontsize=8,va='center')
    b.plot([0,.70],[.25,.25],color=DARK,lw=.6)
    for tick in [0,50,100]:
        x=.7*tick/100;b.plot([x,x],[.25,.22],color=DARK,lw=.6);b.text(x,.18,str(tick),fontsize=7,ha='center')
    b.text(.35,.01,'Complete cases (%)',fontsize=7,ha='center')
    c=panel(fig,[.585,.11,.37,.235],'c','Mouse screen sensitivity',1.16);blank(c)
    ax=fig.add_axes([.67,.215,.225,.085]);clean(ax,grid='x')
    pr=mouse.loc['PRIMARY'];fr=mouse.loc['FITNESS_SENSITIVITY'];cr=mouse.loc['CLASS_CONDITIONED']
    ax.errorbar(float(pr.odds_ratio),1,xerr=[[float(pr.odds_ratio)-float(pr.ci_low)],[float(pr.ci_high)-float(pr.odds_ratio)]],fmt='o',mfc='white',mec=DARK,ecolor=DARK,ms=4,lw=.8,capsize=2)
    ax.scatter(float(fr.fitness_odds_ratio),0,s=20,facecolor='white',edgecolor=DARK,lw=.8)
    ax.axvline(1,color=GRAY,lw=.6,ls=(0,(3,3)));ax.set_yticks([1,0],['Primary','Fitness'],fontsize=7)
    ax.set_ylim(-.5,1.5);ax.set_xlim(0,2.4);ax.set_xticks([0,1,2]);ax.set_xlabel('Odds ratio',fontsize=7)
    c.text(.0,.96,'Primary P = '+pformat(pr.one_sided_fisher_p,3)+'   ·   Fitness P = '+pformat(fr.fitness_one_sided_p,3),fontsize=7)
    ax=fig.add_axes([.67,.09,.225,.046]);clean(ax,grid='x')
    ax.scatter(float(cr.observed_expected_ratio),0,s=20,facecolor='white',edgecolor=DARK,lw=.8)
    ax.axvline(1,color=GRAY,lw=.6,ls=(0,(3,3)));ax.set_xlim(0,2.4);ax.set_xticks([0,1,2]);ax.set_yticks([0],['Class'],fontsize=7)
    ax.set_xlabel('Observed / expected',fontsize=7)
    c.text(0,.27,'Class P = '+pformat(cr.class_conditioned_p,3),fontsize=7)
    save(fig,out,3)

def figure_s4(df,out):
    blood=modules(records(df,'blood_program'));yan=modules(records(df,'yan_program'))
    fig=title(4,'Paired blood scores and Treg-like context',6.65)
    host=panel(fig,[.075,.435,.88,.42],'a','Whole blood: 20 ASD–sibling family pairs',1.08);blank(host)
    # Mean/CI, exact tests and paired directions share row centers.
    ax=fig.add_axes([.16,.46,.33,.34]);clean(ax,grid='x')
    y=np.arange(10);mean=num(blood.mean_difference).to_numpy();lo=num(blood.ci95_low).to_numpy();hi=num(blood.ci95_high).to_numpy()
    ax.errorbar(mean,y,xerr=[mean-lo,hi-mean],fmt='o',mfc='white',mec=NAVY,ecolor=GRAY,ms=4,lw=.9,capsize=2)
    ax.axvline(0,color=GRAY,lw=.7,ls=(0,(3,3)));ax.set_ylim(9.6,-.6);ax.set_yticks(y,['M'+str(i) for i in blood._m],fontsize=8)
    extent=max(abs(lo.min()),abs(hi.max()))*1.13;ax.set_xlim(-extent,extent)
    ax.set_xlabel('Mean ASD − sibling score\n(paired-t 95% display interval)',fontsize=7)
    ax.set_title('Score difference',fontsize=8,pad=13)
    table=fig.add_axes([.535,.46,.20,.34]);table.set_xlim(0,1);table.set_ylim(9.6,-.6);table.axis('off')
    table.text(.2,1.075,'Exact P',transform=table.transAxes,ha='center',fontsize=8)
    table.text(.78,1.075,'BH10 q',transform=table.transAxes,ha='center',fontsize=8)
    directions=fig.add_axes([.78,.46,.17,.34]);directions.set_xlim(0,20);directions.set_ylim(9.6,-.6)
    for j,(_,r) in enumerate(blood.iterrows()):
        table.text(.2,j,pformat(r.p_exact_two_sided,3),ha='center',va='center',fontsize=7.5)
        table.text(.78,j,pformat(r.q_BH_10,3),ha='center',va='center',fontsize=7.5)
        positive,negative,zero=[intval(r[x]) for x in ['n_positive','n_negative','n_zero']]
        assert positive+negative+zero==20
        directions.barh(j,positive,color=NAVY,height=.43,lw=0)
        directions.barh(j,negative,left=positive,color=GRAY,height=.43,lw=0)
        if zero:directions.barh(j,zero,left=positive+negative,color='white',edgecolor=DARK,height=.43,lw=.5)
        directions.text(positive/2,j,str(positive),color='white',ha='center',va='center',fontsize=6.8)
        directions.text(positive+negative/2,j,str(negative),color=DARK,ha='center',va='center',fontsize=6.8)
    directions.set_title('Paired directions',fontsize=8,pad=13);directions.set_yticks([]);directions.set_xticks([0,10,20]);clean(directions,grid=None)
    directions.set_xlabel('')
    for yy,label,col in [(.407,'ASD > sibling',NAVY),(.384,'ASD < sibling',GRAY)]:
        fig.text(.77,yy,'■',fontsize=7,color=col)
        fig.text(.791,yy,label,fontsize=7,color=DARK)
    host=panel(fig,[.075,.09,.88,.195],'b','Yan reported-list overlap',1.30);blank(host)
    for j,subset in enumerate([yan.iloc[:5],yan.iloc[5:]]):
        ax=fig.add_axes([.16+j*.44,.09,.28,.18]);clean(ax,grid='x')
        for k,(_,r) in enumerate(subset.iterrows()):
            value=100*float(r.overlap_fraction)
            ax.barh(k,value,color=GRAY,height=.42,lw=0)
            ax.text(10.8,k,f'{intval(r.program_genes_listed)}/{intval(r.program_size)}',fontsize=7.5,va='center')
        ax.set_yticks(range(5),['M'+str(i) for i in subset._m],fontsize=8);ax.set_ylim(4.6,-.6)
        ax.set_xlim(0,14);ax.set_xticks([0,5,10]);ax.set_xlabel('Program genes in reported list (%)',fontsize=7)
    save(fig,out,4)

def figure_s5(df,out):
    program=modules(records(df,'program'));overlap=records(df,'exact_overlap').set_index('dataset')
    m=overlap.loc['Moreno2026'];y=overlap.loc['Yan2025']
    fig=title(5,'Program support and exact-hit overlap',6.6)
    a=panel(fig,[.075,.435,.875,.41],'a','Moreno formal inference across ten programs',1.08)
    blank(a)
    a.text(0,1.0,'Fixed family: 10 programs   ·   Formal inference: 5 eligible programs',fontsize=8,color=DARK)
    a.text(.02,.89,'Program',fontsize=8,weight='bold');a.text(.22,.89,'Eligibility',fontsize=8,weight='bold')
    a.text(.54,.89,'BH10 q',fontsize=8,weight='bold');a.text(.77,.89,'Formal support',fontsize=8,weight='bold')
    for k,(_,r) in enumerate(program.iterrows()):
        yy=.81-k*.074;eligible=r.formal_eligibility=='ELIGIBLE';supported=r.interpretation_status=='FORMAL_SUPPORTED_Q_LT_005'
        col=TEAL if supported else NAVY if eligible else '#87939D'
        a.text(.02,yy,'M'+str(r._m),fontsize=8,va='center',weight='bold')
        a.text(.22,yy,'Eligible' if eligible else 'Ineligible',fontsize=7.5,va='center',color=col)
        a.text(.54,yy,pformat(r.formal_q_BH10,3) if eligible else '—',fontsize=7.5,va='center',color=col)
        a.scatter(.745,yy,s=18,marker='o' if eligible else 's',facecolor=col if supported or not eligible else 'white',edgecolor=col,lw=.7)
        a.text(.78,yy,'Supported' if supported else 'Not supported' if eligible else 'Not tested',fontsize=7.5,va='center',color=col)
        table_rule(a,yy-.037)
    b=panel(fig,[.085,.09,.37,.225],'b','Exact ASD14 DEG overlap',1.13)
    other_n=intval(m.other_deg)+intval(m.other_non_deg)
    fraction_bars(b,[('ASD14',intval(m.asd14_deg),intval(m.asd14_in_universe)),('Other',intval(m.other_deg),other_n)],6,'DEGs (%)',TEAL)
    b.text(0,.09,f'Fisher enrichment OR = {m.odds_ratio}; P = {m.p_value}',fontsize=7.5)
    b.text(0,-.04,f'{intval(m.tested_universe):,} finite-effect feature rows',fontsize=7.5,color=DARK)
    c=panel(fig,[.575,.10,.375,.21],'c','Measurement and overlap',1.20);blank(c)
    c.text(.28,.92,'Measured',fontsize=7.5,ha='center',weight='bold');c.text(.76,.92,'DEG / list',fontsize=7.5,ha='center',weight='bold')
    c.text(.28,.82,'coverage',fontsize=7.5,ha='center');c.text(.76,.82,'overlap',fontsize=7.5,ha='center')
    table_rule(c,.74)
    for yy,label,coverage,over in [(.58,'Moreno',f'{intval(m.asd14_in_universe)}/{intval(m.asd14_size)}',f'{intval(m.asd14_deg)}/{intval(m.asd14_size)} DEGs'),
                                 (.29,'Yan','Unknown',f'{intval(y.asd14_deg)}/{intval(y.asd14_size)} queried')]:
        c.text(0,yy,label,fontsize=8,weight='bold')
        c.text(.38,yy,coverage,fontsize=7.5,ha='center')
        c.text(.82,yy,over,fontsize=7.5,ha='center')
        table_rule(c,yy-.10)
    save(fig,out,5)

def write_ledger(df,n,out):
    """Every original scalar used by figure families, with source row references."""
    specs={1:{'common_asd_gene':['gene','chen_measured','chen_hit','chen_author_lfc','umhoefer_treg_measured','umhoefer_treg_hit','umhoefer_treg_author_lfc','umhoefer_tconv_measured','umhoefer_tconv_hit','umhoefer_tconv_author_lfc'],
              'tconv_enrichment':['ASD_n','ASD_hits_or_observed','comparator_n','comparator_hits_or_expected','P']},
           2:{'recalibration':['screen','analysis_type_or_removal','A','K','effect','one_sided_P','observed_expected_ratio'],
              'omission':['screen','removed_gene','removed_was_author_hit','one_sided_fisher_p']},
           3:{'class_stratum':['screen','class_label','library_member','N','H','A','K'],'joint_coverage':['screen','COMPLETE_CASE_ASD_MEASURED','ORIGINAL_MEASURED_ASD','COMPLETE_CASE_ASD_COVERAGE'],
              'mouse_sensitivity':['row_type','odds_ratio','ci_low','ci_high','one_sided_fisher_p','fitness_odds_ratio','fitness_one_sided_p','observed_expected_ratio','class_conditioned_p']},
           4:{'blood_program':['program_id','n_pairs','mean_difference','ci95_low','ci95_high','p_exact_two_sided','q_BH_10','n_positive','n_negative','n_zero'],
              'yan_program':['program_id','program_genes_listed','program_size','overlap_fraction']},
           5:{'program':['program_id','formal_eligibility','formal_q_BH10','interpretation_status'],
              'exact_overlap':['dataset','tested_universe','asd14_size','asd14_in_universe','asd14_deg','other_deg','other_non_deg','odds_ratio','p_value']}}
    ledger=[]
    for typ,cols in specs[n].items():
        for _,r in records(df,typ).iterrows():
            for field in cols:
                if r[field]!='':ledger.append({'figure':f'S{n}','record_type':typ,'source_file':r.source_file,'source_row_index':r.source_row_index,'field':field,'exact_value':r[field]})
    qc=out.parent/'04_qc';qc.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(ledger).to_csv(qc/f'S{n}_numerical_value_ledger.tsv',sep='\t',index=False)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--figures',default='1,2,3,4,5');args=ap.parse_args();configure()
    for i in map(int,args.figures.split(',')):
        df=pd.read_csv(args.input_dir/f'Supplementary_Figure_S{i}_source_data.tsv',sep='\t',dtype=str,keep_default_na=False)
        write_ledger(df,i,args.output_dir);globals()[f'figure_s{i}'](df,args.output_dir)
        print(f'RENDERED=Supplementary_Figure_S{i}_FINAL.pdf,png,svg')
if __name__=='__main__':main()
