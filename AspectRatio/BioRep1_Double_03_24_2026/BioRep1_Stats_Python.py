#!/usr/bin/env python3
"""
TIBD Statistical Analysis Suite — Publication-Ready Figures
============================================================
Each analysis produces one clean figure, displayed individually
in PyCharm SciView (or any matplotlib backend).

All values: 2 significant digits.
Style: journal-ready white background, no decorative chrome.

Requirements:
    pip install numpy pandas scipy statsmodels matplotlib seaborn

Usage:
    python tibd_publication.py
    python tibd_publication.py path/to/results.csv
"""

import os, sys, warnings
warnings.filterwarnings('ignore')

import numpy  as np
import pandas as pd
from scipy        import stats
from scipy.stats  import probplot

import statsmodels.formula.api    as smf
from statsmodels.stats.anova      import anova_lm
from statsmodels.stats.multicomp  import pairwise_tukeyhsd
from statsmodels.stats.diagnostic import lilliefors as sm_lilliefors

import matplotlib.pyplot   as plt
import matplotlib.patches  as mpatches
import matplotlib.ticker   as mticker
import matplotlib.gridspec as gridspec
import seaborn as sns

# ═══════════════════════════════════════════════════════════════
#  CONFIGURATION
# ═══════════════════════════════════════════════════════════════
CSV_PATH   = (
    '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio'
    '/BioRep1_Double_03_24_2026/BioRep1_Full_Unfiltered_Results.csv'
)
COND_ORDER = ['TC', 'Random', 'Aligned']

# ═══════════════════════════════════════════════════════════════
#  PUBLICATION STYLE
# ═══════════════════════════════════════════════════════════════
# Two-color palette — blue & orange (colorblind-safe)
BLUE   = '#2166AC'
ORANGE = '#D6604D'
PAL    = [BLUE, ORANGE]

GREY_DARK  = '#222222'
GREY_MID   = '#555555'
GREY_LIGHT = '#AAAAAA'
GREY_GRID  = '#E5E5E5'

plt.rcParams.update({
    # Canvas
    'figure.facecolor'     : 'white',
    'figure.dpi'           : 150,
    # Axes
    'axes.facecolor'       : 'white',
    'axes.edgecolor'       : GREY_DARK,
    'axes.linewidth'       : 0.8,
    'axes.labelcolor'      : GREY_DARK,
    'axes.labelsize'       : 8,
    'axes.titlesize'       : 9,
    'axes.titlepad'        : 6,
    'axes.titleweight'     : 'bold',
    'axes.titlecolor'      : GREY_DARK,
    'axes.spines.top'      : False,
    'axes.spines.right'    : False,
    # Ticks
    'xtick.color'          : GREY_DARK,
    'ytick.color'          : GREY_DARK,
    'xtick.labelsize'      : 7.5,
    'ytick.labelsize'      : 7.5,
    'xtick.major.width'    : 0.7,
    'ytick.major.width'    : 0.7,
    'xtick.major.size'     : 3,
    'ytick.major.size'     : 3,
    # Grid
    'grid.color'           : GREY_GRID,
    'grid.linewidth'       : 0.6,
    # Legend
    'legend.facecolor'     : 'white',
    'legend.edgecolor'     : GREY_LIGHT,
    'legend.fontsize'      : 7.5,
    'legend.title_fontsize': 8,
    'legend.framealpha'    : 0.9,
    # Lines & markers
    'lines.linewidth'      : 1.5,
    'lines.antialiased'    : True,
    'patch.antialiased'    : True,
    # Text
    'text.color'           : GREY_DARK,
    'font.family'          : ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size'            : 8,
    # Save
    'savefig.dpi'          : 300,
    'savefig.bbox'         : 'tight',
    'savefig.facecolor'    : 'white',
})

# ═══════════════════════════════════════════════════════════════
#  FORMATTING HELPERS
# ═══════════════════════════════════════════════════════════════
def sf(x, n=2):
    """Format a number to n significant figures."""
    if pd.isna(x) or not np.isfinite(float(x)):
        return 'n/a'
    if x == 0:
        return '0'
    mag = int(np.floor(np.log10(abs(x))))
    rounded = round(x, -mag + (n - 1))
    if abs(mag) >= 3 or mag < -1:
        return f'{x:.{n-1}e}'
    decimals = max(0, n - 1 - mag)
    return f'{rounded:.{decimals}f}'

def pval(p, n=2):
    """Format a p-value to n significant figures, e.g. 2.3e-04."""
    if pd.isna(p): return 'n/a'
    if p < 0.0001:
        return f'{p:.{n-1}e}'
    return f'{p:.{n}g}'

def stars(p):
    if pd.isna(p): return ''
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

def clean_ax(ax, grid_axis='y'):
    """Minimal publication spine + grid."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(GREY_DARK)
    ax.spines['bottom'].set_color(GREY_DARK)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GREY_GRID, linewidth=0.6, zorder=0)
    ax.tick_params(colors=GREY_DARK, length=3, width=0.7)

def sig_bracket(ax, x1, x2, y, p, h=0.015, dy=0.005):
    """Draw a significance bracket above two bars/groups."""
    yy = y + h
    ax.plot([x1, x1, x2, x2], [y, yy, yy, y],
            lw=0.8, color=GREY_DARK, clip_on=False)
    lbl = stars(p) if stars(p) != 'ns' else 'ns'
    ax.text((x1+x2)/2, yy + dy, lbl,
            ha='center', va='bottom', fontsize=7,
            color=GREY_DARK)

# ═══════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════
def load_data(path):
    df = pd.read_csv(path)
    df['CellLine']  = df['CellLine'].astype(str)
    df['Condition'] = df['Condition'].astype(str)
    df = df.dropna(subset=['AspectRatio'])
    df['Condition'] = pd.Categorical(df['Condition'],
                                     categories=COND_ORDER, ordered=True)
    df['Sample'] = df['Sample'].astype(str)
    df['group']  = df['CellLine'] + ' | ' + df['Condition'].astype(str)
    return df

# ═══════════════════════════════════════════════════════════════
#  STATISTICS
# ═══════════════════════════════════════════════════════════════
def run_stats(df):
    R = {}
    R['cl']   = sorted(df['CellLine'].unique())
    R['cond'] = list(COND_ORDER)

    # Two-Way ANOVA
    model = smf.ols(
        'Q("AspectRatio") ~ Q("CellLine") * Q("Condition")', data=df
    ).fit()
    atbl = anova_lm(model, typ=2)
    R['model'] = model
    R['atbl']  = atbl

    k_cl  = 'Q("CellLine")';  k_cd = 'Q("Condition")'
    k_int = 'Q("CellLine"):Q("Condition")';  k_res = 'Residual'

    ss = {k: atbl.loc[v, 'sum_sq']
          for k, v in [('cl',k_cl),('cond',k_cd),('int',k_int),('err',k_res)]}
    ss['tot'] = sum(ss.values())
    dfd = {k: atbl.loc[v, 'df']
           for k, v in [('cl',k_cl),('cond',k_cd),('int',k_int),('err',k_res)]}
    Fv  = {k: atbl.loc[v, 'F']       for k,v in [('cl',k_cl),('cond',k_cd),('int',k_int)]}
    pv  = {k: atbl.loc[v, 'PR(>F)']  for k,v in [('cl',k_cl),('cond',k_cd),('int',k_int)]}

    R['ss']=ss; R['df']=dfd; R['F']=Fv; R['p']=pv
    R['eta2']  = {k: ss[k]/ss['tot']          for k in ['cl','cond','int']}
    R['peta2'] = {k: ss[k]/(ss[k]+ss['err'])  for k in ['cl','cond','int']}

    # Normality (Lilliefors per group)
    norm = []
    for cl in R['cl']:
        for cd in R['cond']:
            d = df.loc[(df['CellLine']==cl)&(df['Condition']==cd),'AspectRatio'].values
            stat, p = sm_lilliefors(d) if len(d) >= 4 else (np.nan, np.nan)
            norm.append({'label': f'{cl} | {cd}', 'n': len(d),
                         'stat': stat, 'p': p, 'data': d, 'cl': cl, 'cond': cd})
    R['norm'] = norm

    # Levene's
    active = sorted([g for g in df['group'].unique() if pd.notna(g)])
    grps   = [df.loc[df['group']==g,'AspectRatio'].values for g in active]
    R['levene']    = dict(zip(['stat','p'], stats.levene(*grps, center='median') if len(grps)>1 else (np.nan,np.nan)))
    R['grp_stats'] = [{'group':g,'n':len(d),'mean':d.mean(),'std':d.std(),'var':d.var()}
                      for g,d in zip(active,grps)]

    # Means & SEMs
    nc, nd = len(R['cl']), len(R['cond'])
    means = np.zeros((nc,nd)); sems = np.zeros((nc,nd))
    for i,cl in enumerate(R['cl']):
        for j,cd in enumerate(R['cond']):
            d = df.loc[(df['CellLine']==cl)&(df['Condition']==cd),'AspectRatio'].values
            means[i,j]=d.mean(); sems[i,j]=d.std()/np.sqrt(len(d))
    R['means']=means; R['sems']=sems

    # Tukey-Kramer
    mask  = df['group'].notna() & (df['group'].astype(str)!='nan')
    cdf   = df[mask]
    R['tukey'] = pairwise_tukeyhsd(cdf['AspectRatio'], cdf['group'].astype(str))

    # Linear Mixed-Effects
    try:
        R['lme'] = smf.mixedlm(
            'AspectRatio ~ C(CellLine) * C(Condition)',
            data=df, groups=df['Sample']
        ).fit(reml=True)
    except Exception as e:
        R['lme'] = None; R['lme_err'] = str(e)

    return R


# ═══════════════════════════════════════════════════════════════
#  FIGURE 1 — TWO-WAY ANOVA TABLE
# ═══════════════════════════════════════════════════════════════
def fig_anova(R):
    fig, ax = plt.subplots(figsize=(7, 2.8))
    ax.axis('off')

    ms = lambda s, d: s/d
    total_df = int(round(sum(R['df'].values())))
    col_labels = ['Source', 'SS', 'df', 'MS', 'F', 'p', '']
    col_w      = [0.26, 0.12, 0.06, 0.12, 0.10, 0.14, 0.10]

    rows = [
        ['Cell Line',
         sf(R['ss']['cl']),    str(int(R['df']['cl'])),
         sf(ms(R['ss']['cl'],   R['df']['cl'])),
         sf(R['F']['cl']),     pval(R['p']['cl']), stars(R['p']['cl'])],
        ['Condition',
         sf(R['ss']['cond']),  str(int(R['df']['cond'])),
         sf(ms(R['ss']['cond'],R['df']['cond'])),
         sf(R['F']['cond']),   pval(R['p']['cond']), stars(R['p']['cond'])],
        ['Cell Line × Condition',
         sf(R['ss']['int']),   str(int(R['df']['int'])),
         sf(ms(R['ss']['int'], R['df']['int'])),
         sf(R['F']['int']),    pval(R['p']['int']), stars(R['p']['int'])],
        ['Error',
         sf(R['ss']['err']),   str(int(R['df']['err'])),
         sf(ms(R['ss']['err'],R['df']['err'])),
         '—', '—', '—'],
        ['Total',
         sf(R['ss']['tot']),   str(total_df),
         '—', '—', '—', '—'],
    ]

    t = ax.table(
        cellText=rows, colLabels=col_labels,
        colWidths=col_w,
        cellLoc='center', loc='center',
    )
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1, 1.9)

    # Header row
    for j in range(len(col_labels)):
        cell = t[0, j]
        cell.set_facecolor('#DDEEFF')
        cell.set_text_props(fontweight='bold', color=GREY_DARK)
        cell.set_edgecolor('white')
    # Data rows
    for i in range(len(rows)):
        row_bg = 'white' if i%2==0 else '#F7F9FC'
        sig_row = (i < 3 and R['p'][['cl','cond','int'][i]] < 0.05)
        for j in range(len(col_labels)):
            cell = t[i+1, j]
            cell.set_facecolor('#FFF5F5' if sig_row else row_bg)
            cell.set_edgecolor('white')
            cell.set_text_props(color=GREY_DARK)

    ax.set_title('Two-Way ANOVA — Aspect Ratio', pad=10, loc='left',
                 fontsize=9, fontweight='bold')
    fig.text(0.02, 0.02, '* p < 0.05   ** p < 0.01   *** p < 0.001',
             fontsize=6.5, color=GREY_MID)
    fig.tight_layout()
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 2 — NORMALITY (QQ PLOTS)
# ═══════════════════════════════════════════════════════════════
def fig_normality(R):
    norm = R['norm']
    nCL  = len(R['cl']); nCd = len(R['cond'])

    fig, axes = plt.subplots(nCL, nCd,
                              figsize=(2.5*nCd, 2.5*nCL),
                              squeeze=False)
    fig.suptitle('Normality — Q–Q Plots per Group',
                 fontsize=9, fontweight='bold', y=1.01, x=0.04, ha='left')

    idx = 0
    for i in range(nCL):
        for j in range(nCd):
            ax  = axes[i][j]
            nd  = norm[idx]
            if len(nd['data']) >= 4:
                (osm, osr), (slope, intercept, _) = probplot(nd['data'], dist='norm')
                ax.scatter(osm, osr, s=1.5, color=PAL[i],
                           alpha=0.3, rasterized=True, linewidths=0)
                ax.plot(osm, slope*np.array(osm)+intercept,
                        color=GREY_DARK, lw=1.2, alpha=0.85)
            clean_ax(ax, grid_axis=None)
            title = f'{nd["cl"]} | {nd["cond"]}\nn = {nd["n"]:,}'
            ax.set_title(title, fontsize=7.5, pad=4, loc='left')
            ax.set_xlabel('Theoretical quantiles', fontsize=7)
            ax.set_ylabel('Sample quantiles',       fontsize=7)
            if not pd.isna(nd['p']):
                vc  = '#CC0000' if nd['p'] < 0.05 else '#006622'
                lbl = f'p = {pval(nd["p"])}'
                ax.text(0.97, 0.04, lbl, transform=ax.transAxes,
                        ha='right', va='bottom', fontsize=7,
                        color=vc, fontweight='bold')
            idx += 1

    fig.tight_layout(rect=[0,0,1,0.98])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 3 — LEVENE'S TEST
# ═══════════════════════════════════════════════════════════════
def fig_levene(R):
    lev = R['levene']
    gs  = R['grp_stats']

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5),
                              gridspec_kw={'width_ratios': [1, 1.4]})
    fig.suptitle("Levene's Test — Homogeneity of Variance",
                 fontsize=9, fontweight='bold', x=0.02, ha='left')

    # ── Variance bar chart ──────────────────────────────────────
    ax = axes[0]
    labels = [g['group'].split('|')[-1].strip() for g in gs]
    cl_ids = [g['group'].split('|')[0].strip() for g in gs]
    vars_  = [g['var'] for g in gs]
    cols   = [PAL[0] if c == R['cl'][0] else PAL[1] for c in cl_ids]
    ypos   = range(len(labels))
    bars   = ax.barh(list(ypos), vars_, color=cols, alpha=0.82,
                     edgecolor='white', linewidth=0.8, height=0.55)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel('Variance', fontsize=8)
    for bar, val in zip(bars, vars_):
        ax.text(val + max(vars_)*0.02, bar.get_y()+bar.get_height()/2,
                sf(val), va='center', fontsize=7, color=GREY_MID)
    clean_ax(ax, grid_axis='x')
    ax.grid(True, axis='x', color=GREY_GRID, linewidth=0.6)
    ax.grid(False, axis='y')
    ax.set_title('Variance per Group', fontsize=8.5, fontweight='bold', loc='left')

    lp = lev['p']
    vc = '#CC0000' if lp < 0.05 else '#006622'
    vt = 'Heterogeneous' if lp < 0.05 else 'Homogeneous'
    ax.text(0.98, 0.02, f'{vt}  (p = {pval(lp)})',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=7.5, color=vc, fontweight='bold')

    # ── Descriptive table ───────────────────────────────────────
    ax2 = axes[1]
    ax2.axis('off')
    col_labels = ['Group', 'n', 'Mean', 'SD', 'Variance']
    col_w      = [0.32, 0.12, 0.16, 0.16, 0.22]
    tbl_rows   = [
        [g['group'], f"{g['n']:,}", sf(g['mean']), sf(g['std']), sf(g['var'])]
        for g in gs
    ]
    t = ax2.table(cellText=tbl_rows, colLabels=col_labels,
                  colWidths=col_w, cellLoc='center', loc='center')
    t.auto_set_font_size(False); t.set_fontsize(8); t.scale(1, 1.8)
    for j in range(len(col_labels)):
        c = t[0,j]; c.set_facecolor('#DDEEFF')
        c.set_text_props(fontweight='bold'); c.set_edgecolor('white')
    for i in range(len(tbl_rows)):
        bg = 'white' if i%2==0 else '#F7F9FC'
        for j in range(len(col_labels)):
            c = t[i+1,j]; c.set_facecolor(bg); c.set_edgecolor('white')
    ax2.set_title('Descriptive Statistics', fontsize=8.5,
                  fontweight='bold', loc='left', pad=4)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 4 — INTERACTION PLOT
# ═══════════════════════════════════════════════════════════════
def fig_interaction(R):
    x = np.arange(len(R['cond']))

    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    for i, cl in enumerate(R['cl']):
        m, s = R['means'][i], R['sems'][i]
        ax.fill_between(x, m-s, m+s, color=PAL[i], alpha=0.10)
        ax.plot(x, m, '-o', color=PAL[i], label=cl,
                lw=1.8, ms=7,
                markerfacecolor=PAL[i], markeredgecolor='white', markeredgewidth=1.4)
        ax.errorbar(x, m, yerr=s, fmt='none', ecolor=PAL[i],
                    elinewidth=1.2, capsize=4, capthick=1.2, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(R['cond'], fontsize=8)
    ax.set_xlabel('Scaffold Condition', fontsize=8)
    ax.set_ylabel('Mean Aspect Ratio ± SEM', fontsize=8)
    ax.legend(title='Cell Line', framealpha=0.9,
              loc='best', fontsize=7.5, title_fontsize=8)
    clean_ax(ax)

    p_int = R['p']['int']
    ax.set_title(f'Interaction: Cell Line × Condition\n'
                 f'p = {pval(p_int)}  {stars(p_int)}',
                 fontsize=9, fontweight='bold', loc='left')
    fig.tight_layout()
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 5 — TUKEY-KRAMER POST-HOC
# ═══════════════════════════════════════════════════════════════
def fig_posthoc(R):
    tukey = R['tukey']
    raw   = tukey.summary().data[1:]
    g1    = [str(x[0]) for x in raw]
    g2    = [str(x[1]) for x in raw]
    md    = [float(x[2]) for x in raw]
    padj  = [float(x[3]) for x in raw]
    lo    = [float(x[4]) for x in raw]
    hi    = [float(x[5]) for x in raw]
    rej   = [str(x[6]).lower()=='true' for x in raw]
    n     = len(g1)

    fig, axes = plt.subplots(1, 2, figsize=(14, max(4, n*0.4+1.2)),
                              gridspec_kw={'width_ratios': [1, 1.5]})
    fig.suptitle('Post-Hoc — Tukey–Kramer Pairwise Comparisons',
                 fontsize=9, fontweight='bold', x=0.02, ha='left')

    # ── Forest plot ─────────────────────────────────────────────
    ax = axes[0]
    y_pos = np.arange(n, 0, -1)
    for i in range(n):
        col   = BLUE if rej[i] else GREY_LIGHT
        alpha = 1.0  if rej[i] else 0.6
        ax.plot([lo[i], hi[i]], [y_pos[i]]*2, '-', color=col, lw=2.0, alpha=alpha)
        ax.plot(md[i], y_pos[i], 'o', color=col, ms=5,
                markeredgecolor='white', mew=1.0, alpha=alpha)
    ax.axvline(0, color=GREY_DARK, lw=0.9, ls='--', alpha=0.6)
    short = [f"{g1[i]} vs {g2[i]}" for i in range(n)]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(short, fontsize=7.5)
    ax.set_xlabel('Mean Difference (95% CI)', fontsize=8)
    ax.set_title('Comparison Intervals', fontsize=8.5, fontweight='bold', loc='left')
    clean_ax(ax, grid_axis='x')
    ax.grid(True, axis='x', color=GREY_GRID, linewidth=0.6)
    ax.grid(False, axis='y')
    ax.legend(handles=[
        mpatches.Patch(color=BLUE,      label='p < 0.05'),
        mpatches.Patch(color=GREY_LIGHT,label='Not significant'),
    ], fontsize=7.5, loc='lower right')

    # ── Results table ────────────────────────────────────────────
    ax2 = axes[1]
    ax2.axis('off')
    col_labels = ['Group 1', 'Group 2', 'Δ Mean', 'Adj. p', '']
    tbl_rows = [
        [g1[i],
         g2[i],
         sf(md[i]), pval(padj[i]), stars(padj[i])]
        for i in range(n)
    ]
    t = ax2.table(cellText=tbl_rows, colLabels=col_labels,
                  colWidths=[0.35,0.35,0.10,0.10,0.10],
                  cellLoc='center', loc='center')
    t.auto_set_font_size(False); t.set_fontsize(8); t.scale(1, 1.65)
    for j in range(len(col_labels)):
        c = t[0,j]; c.set_facecolor('#DDEEFF')
        c.set_text_props(fontweight='bold'); c.set_edgecolor('white')
    for i, row in enumerate(tbl_rows):
        bg = 'white' if i%2==0 else '#F7F9FC'
        sig_row = rej[i]
        for j in range(len(col_labels)):
            c = t[i+1,j]
            c.set_facecolor('#EFF5FF' if sig_row else bg)
            c.set_edgecolor('white')
            c.set_text_props(color=BLUE if sig_row and j==4 else GREY_DARK,
                             fontweight='bold' if sig_row and j==4 else 'normal')
    ax2.set_title('All Pairwise Results', fontsize=8.5,
                  fontweight='bold', loc='left', pad=4)

    fig.tight_layout(rect=[0,0,1,0.94])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 6 — EFFECT SIZES
# ═══════════════════════════════════════════════════════════════
def fig_effects(R):
    sources = ['Cell Line', 'Condition', 'CL × Cond']
    ev   = [R['eta2']['cl'],  R['eta2']['cond'],  R['eta2']['int']]
    pev  = [R['peta2']['cl'], R['peta2']['cond'], R['peta2']['int']]

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.0),
                              gridspec_kw={'width_ratios':[1.6, 1]})
    fig.suptitle('Effect Sizes — η² and Partial η²',
                 fontsize=9, fontweight='bold', x=0.02, ha='left')

    # ── Bar chart ────────────────────────────────────────────────
    ax = axes[0]
    x = np.arange(3); w = 0.32
    b1 = ax.bar(x-w/2, ev,  w, label='η²', color=BLUE,   alpha=0.82,
                edgecolor='white', lw=0.8, zorder=4)
    b2 = ax.bar(x+w/2, pev, w, label='Partial η²', color=ORANGE, alpha=0.82,
                edgecolor='white', lw=0.8, zorder=4)
    for bar, val in zip(list(b1)+list(b2), ev+pev):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.001,
                sf(val), ha='center', va='bottom', fontsize=7.5)
    for thr, lbl, ls in [(0.01,'Small',':'),(0.06,'Medium','--'),(0.15,'Large','-')]:
        ax.axhline(thr, color=GREY_LIGHT, lw=0.9, ls=ls)
        ax.text(2.66, thr+0.002, lbl, fontsize=7, color=GREY_MID)
    ax.set_xticks(x)
    ax.set_xticklabels(sources, fontsize=8)
    ax.set_ylabel('Effect Size', fontsize=8)
    ax.legend(fontsize=7.5, loc='upper left', framealpha=0.9)
    clean_ax(ax)

    # ── Summary table ────────────────────────────────────────────
    ax2 = axes[1]; ax2.axis('off')
    def emag(pe):
        return 'Large' if pe>=.14 else 'Medium' if pe>=.06 else 'Small' if pe>=.01 else 'Negligible'
    t = ax2.table(
        cellText=[[s, sf(e), sf(pe), emag(pe)] for s,e,pe in zip(sources,ev,pev)],
        colLabels=['Source','η²','Partial η²','Magnitude'],
        colWidths=[0.30,0.17,0.22,0.31],
        cellLoc='center', loc='center',
    )
    t.auto_set_font_size(False); t.set_fontsize(8); t.scale(1, 2.0)
    for j in range(4):
        c=t[0,j]; c.set_facecolor('#DDEEFF')
        c.set_text_props(fontweight='bold'); c.set_edgecolor('white')
    for i in range(3):
        bg = 'white' if i%2==0 else '#F7F9FC'
        for j in range(4):
            c=t[i+1,j]; c.set_facecolor(bg); c.set_edgecolor('white')
    ax2.set_title('Summary', fontsize=8.5, fontweight='bold', loc='left', pad=4)

    fig.tight_layout(rect=[0,0,1,0.94])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 7 — LINEAR MIXED-EFFECTS MODEL
# ═══════════════════════════════════════════════════════════════
def fig_lme(R):
    lme = R.get('lme')

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8),
                              gridspec_kw={'width_ratios':[1.8, 1]})
    fig.suptitle(
        'Linear Mixed-Effects — AspectRatio ~ CellLine × Condition + (1 | Sample)',
        fontsize=9, fontweight='bold', x=0.02, ha='left',
    )

    if lme is None:
        for ax in axes: ax.axis('off')
        axes[0].text(0.5, 0.5, f"Fit failed: {R.get('lme_err','unknown')}",
                     ha='center', va='center', color='#CC0000', fontsize=9)
        plt.show(); return

    # ── Fixed effects table ──────────────────────────────────────
    ax = axes[0]; ax.axis('off')
    fe_rows = []
    for pn in lme.params.index:
        short = (pn.replace('C(CellLine)','CL:').replace('C(Condition)','Cond:')
                   .replace('[T.','[').replace('Intercept','Intercept'))
        ci = lme.conf_int().loc[pn]
        fe_rows.append([
            short,
            sf(lme.params[pn]), sf(lme.bse[pn]),
            sf(lme.tvalues[pn]), pval(lme.pvalues[pn]),
            stars(lme.pvalues[pn]),
            f'[{sf(ci[0])}, {sf(ci[1])}]',
        ])
    t = ax.table(
        cellText=fe_rows,
        colLabels=['Parameter','Coef.','SE','t','p','','95% CI'],
        colWidths=[0.24,0.10,0.09,0.09,0.12,0.08,0.28],
        cellLoc='center', loc='center',
    )
    t.auto_set_font_size(False); t.set_fontsize(7.8); t.scale(1,1.75)
    for j in range(7):
        c=t[0,j]; c.set_facecolor('#DDEEFF')
        c.set_text_props(fontweight='bold'); c.set_edgecolor('white')
    for i in range(len(fe_rows)):
        bg = 'white' if i%2==0 else '#F7F9FC'
        sig_row = lme.pvalues.iloc[i] < 0.05
        for j in range(7):
            c=t[i+1,j]
            c.set_facecolor('#EFF5FF' if sig_row else bg)
            c.set_edgecolor('white')
    ax.set_title('Fixed Effects', fontsize=8.5, fontweight='bold', loc='left', pad=4)

    # ── Model fit callout ────────────────────────────────────────
    ax2 = axes[1]; ax2.axis('off')
    try:
        aic = lme.aic; bic = lme.bic; llf = lme.llf
    except Exception:
        aic = bic = llf = np.nan
    re_var  = lme.cov_re.values[0][0] if lme.cov_re is not None else np.nan
    res_var = lme.scale
    icc = re_var/(re_var+res_var) if not any(np.isnan([re_var, res_var])) else np.nan

    fit_rows = [
        ['Log-Likelihood', sf(llf)],
        ['AIC',            sf(aic)],
        ['BIC',            sf(bic)],
        ['Sample Var.',    sf(re_var)],
        ['Residual Var.',  sf(res_var)],
        ['ICC',            sf(icc)],
    ]
    t2 = ax2.table(cellText=fit_rows,
                   colLabels=['Fit Index','Value'],
                   colWidths=[0.58, 0.42],
                   cellLoc='center', loc='center')
    t2.auto_set_font_size(False); t2.set_fontsize(8); t2.scale(1, 1.9)
    for j in range(2):
        c=t2[0,j]; c.set_facecolor('#DDEEFF')
        c.set_text_props(fontweight='bold'); c.set_edgecolor('white')
    for i in range(len(fit_rows)):
        bg = 'white' if i%2==0 else '#F7F9FC'
        for j in range(2):
            c=t2[i+1,j]; c.set_facecolor(bg); c.set_edgecolor('white')
    ax2.set_title('Model Fit', fontsize=8.5, fontweight='bold', loc='left', pad=4)

    fig.tight_layout(rect=[0,0,1,0.94])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 8 — BOX PLOTS
# ═══════════════════════════════════════════════════════════════
def fig_boxplots(df, R):
    nCL = len(R['cl'])
    fig, axes = plt.subplots(1, nCL, figsize=(3.2*nCL, 3.8), squeeze=False)
    fig.suptitle('Box Plots — Aspect Ratio by Condition',
                 fontsize=9, fontweight='bold', x=0.02, ha='left')

    for i, cl in enumerate(R['cl']):
        ax  = axes[0][i]
        sub = df[df['CellLine']==cl]
        cond_data = [sub[sub['Condition']==cd]['AspectRatio'].values
                     for cd in R['cond']]
        bp = ax.boxplot(
            cond_data,
            patch_artist=True, widths=0.44,
            medianprops  =dict(color='white', lw=2.0),
            whiskerprops =dict(color=GREY_MID, lw=1.0),
            capprops     =dict(color=GREY_MID, lw=1.0),
            flierprops   =dict(marker='o', ms=2.5, alpha=0.3,
                               markerfacecolor=PAL[i], markeredgecolor='none'),
            boxprops     =dict(linewidth=0),
            showfliers=False
        )
        for patch in bp['boxes']:
            patch.set_facecolor(PAL[i]); patch.set_alpha(0.78)
        ax.set_xticks(range(1, len(R['cond'])+1))
        ax.set_xticklabels(R['cond'], fontsize=8)
        ax.set_title(cl, fontsize=9, fontweight='bold', loc='left')
        ax.set_xlabel('Condition', fontsize=8)
        ax.set_ylabel('Aspect Ratio', fontsize=8) if i==0 else None
        clean_ax(ax)

    fig.tight_layout(rect=[0,0,1,0.94])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 9 — VIOLIN PLOTS
# ═══════════════════════════════════════════════════════════════
def fig_violins(df, R):
    nCL = len(R['cl'])
    fig, axes = plt.subplots(1, nCL, figsize=(3.2*nCL, 3.8), squeeze=False)
    fig.suptitle('Violin Plots — Distribution Shape by Condition',
                 fontsize=9, fontweight='bold', x=0.02, ha='left')

    for i, cl in enumerate(R['cl']):
        ax  = axes[0][i]
        sub = df[df['CellLine']==cl]
        viol_data = [sub[sub['Condition']==cd]['AspectRatio'].values
                     for cd in R['cond']]
        parts = ax.violinplot(viol_data,
                              positions=range(1, len(R['cond'])+1),
                              showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(PAL[i])
            pc.set_edgecolor(GREY_MID)
            pc.set_alpha(0.70)
            pc.set_linewidth(0.6)
        parts['cmedians'].set_color(GREY_DARK)
        parts['cmedians'].set_linewidth(1.8)
        ax.set_xticks(range(1, len(R['cond'])+1))
        ax.set_xticklabels(R['cond'], fontsize=8)
        ax.set_title(cl, fontsize=9, fontweight='bold', loc='left')
        ax.set_xlabel('Condition', fontsize=8)
        ax.set_ylabel('Aspect Ratio', fontsize=8) if i==0 else None
        clean_ax(ax)
        ax.set_ylim(top=5)

    fig.tight_layout(rect=[0,0,1,0.94])
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  FIGURE 10 — MEAN ± SEM BAR CHART
# ═══════════════════════════════════════════════════════════════
def fig_meansem(R):
    nCL = len(R['cl'])
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    x = np.arange(len(R['cond'])); w = 0.32
    for i, cl in enumerate(R['cl']):
        offset = (i - (nCL-1)/2) * w
        ax.bar(x + offset, R['means'][i], w,
               label=cl, color=PAL[i], alpha=0.82,
               edgecolor='white', lw=0.8, zorder=4)
        ax.errorbar(x + offset, R['means'][i], yerr=R['sems'][i],
                    fmt='none', ecolor=GREY_DARK,
                    elinewidth=1.2, capsize=4, capthick=1.2, zorder=5)

    ax.set_xticks(x)
    ax.set_xticklabels(R['cond'], fontsize=8)
    ax.set_xlabel('Scaffold Condition', fontsize=8)
    ax.set_ylabel('Mean Aspect Ratio ± SEM', fontsize=8)
    ax.legend(title='Cell Line', fontsize=7.5, title_fontsize=8,
              loc='best', framealpha=0.9)
    clean_ax(ax)
    ax.set_title('Mean ± SEM', fontsize=9, fontweight='bold', loc='left')

    fig.tight_layout()
    plt.show()


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else CSV_PATH
    if not os.path.exists(path):
        print(f'[ERROR] CSV not found:\n  {path}')
        sys.exit(1)

    print(f'Loading: {path}')
    df = load_data(path)
    print(f'  {len(df):,} cells  |  {len(df["group"].unique())} groups')

    print('Running statistics …')
    R = run_stats(df)
    print('  Done.\n')

    print('Rendering figures (each will appear separately in SciView) …')
    fig_anova(R)          # Figure 1
    fig_normality(R)      # Figure 2
    fig_levene(R)         # Figure 3
    fig_interaction(R)    # Figure 4
    fig_posthoc(R)        # Figure 5
    fig_effects(R)        # Figure 6
    fig_lme(R)            # Figure 7
    fig_boxplots(df, R)   # Figure 8
    fig_violins(df, R)    # Figure 9
    fig_meansem(R)        # Figure 10
    print('All figures rendered.')