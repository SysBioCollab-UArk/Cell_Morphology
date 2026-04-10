#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════╗
║         TIBD  Statistical Analysis Suite  —  Python Edition     ║
║  Ten-tab pipeline:  ANOVA · Normality · Levene's · Interaction  ║
║  Post-Hoc · Effect Sizes · Mixed-Effects · Box · Violin · SEM   ║
╚══════════════════════════════════════════════════════════════════╝

Requirements (pip install):
    numpy pandas scipy statsmodels matplotlib seaborn

Usage:
    python tibd_analysis.py                        # uses hard-coded path below
    python tibd_analysis.py path/to/results.csv    # or pass CSV as argument
"""

import os, sys, warnings
warnings.filterwarnings('ignore')

import numpy  as np
import pandas as pd
from scipy        import stats
from scipy.stats  import probplot

import statsmodels.formula.api             as smf
from statsmodels.stats.anova               import anova_lm
from statsmodels.stats.multicomp           import pairwise_tukeyhsd
from statsmodels.stats.diagnostic          import lilliefors as sm_lilliefors

import tkinter as tk
from tkinter import ttk, filedialog

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot   as plt
import matplotlib.patches  as mpatches
import matplotlib.ticker   as mticker
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import seaborn as sns

# ══════════════════════════════════════════════════════════════════
#  CONFIGURATION  ← edit CSV_PATH to match your file location
# ══════════════════════════════════════════════════════════════════
CSV_PATH   = (
    '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio'
    '/BioRep1_Double_03_24_2026/BioRep1_Final_Results.csv'
)
COND_ORDER = ['TC', 'Random', 'Aligned']

# ══════════════════════════════════════════════════════════════════
#  DESIGN TOKENS
# ══════════════════════════════════════════════════════════════════
C = dict(
    # primary palette
    blue   = '#2980B9',
    orange = '#E67E22',
    # UI chrome
    bg     = '#0F1923',      # near-black app background
    panel  = '#18263A',      # card/panel background
    card   = '#1E3050',      # slightly lighter card
    hdr    = '#0A1628',      # darkest — tab bar / header
    # text
    bright = '#E8F4FD',
    text   = '#A8C0D6',
    muted  = '#5E7A99',
    # semantics
    red    = '#E74C3C',
    green  = '#2ECC71',
    amber  = '#F39C12',
    # plot chrome
    grid   = '#1E3050',
    border = '#2C4A6E',
    # table
    th_bg  = '#0D2137',
    row_e  = '#162030',
    row_o  = '#1C2D44',
)
PAL = [C['blue'], C['orange']]

# matplotlib style applied globally to all figures
plt.rcParams.update({
    'figure.facecolor'  : C['panel'],
    'axes.facecolor'    : C['panel'],
    'axes.edgecolor'    : C['border'],
    'axes.labelcolor'   : C['text'],
    'axes.titlecolor'   : C['bright'],
    'xtick.color'       : C['text'],
    'ytick.color'       : C['text'],
    'text.color'        : C['text'],
    'grid.color'        : C['grid'],
    'grid.linewidth'    : 0.6,
    'legend.facecolor'  : C['card'],
    'legend.edgecolor'  : C['border'],
    'legend.labelcolor' : C['text'],
    'font.family'       : 'DejaVu Sans',
    'axes.spines.top'   : False,
    'axes.spines.right' : False,
})

# ══════════════════════════════════════════════════════════════════
#  HELPERS
# ══════════════════════════════════════════════════════════════════
def sig(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

def emag(pe):
    return 'LARGE' if pe >= 0.14 else 'MEDIUM' if pe >= 0.06 else 'SMALL' if pe >= 0.01 else 'NEGLIGIBLE'

def polish(ax, title='', xl='', yl='', fs=12):
    ax.set_facecolor(C['panel'])
    ax.spines['left'].set_color(C['border'])
    ax.spines['bottom'].set_color(C['border'])
    ax.tick_params(colors=C['text'], labelsize=9.5)
    ax.grid(True, color=C['grid'], linewidth=0.6, zorder=0)
    if title: ax.set_title(title, color=C['bright'], fontsize=fs,
                            fontweight='bold', pad=10, loc='left')
    if xl:    ax.set_xlabel(xl, color=C['text'], fontsize=10.5, labelpad=6)
    if yl:    ax.set_ylabel(yl, color=C['text'], fontsize=10.5, labelpad=6)

def draw_table(ax, hdrs, rows, cw=None, title=''):
    """Render a styled dark-theme table inside a matplotlib axes."""
    ax.set_facecolor(C['panel'])
    ax.axis('off')
    if not cw:
        cw = [1.0 / len(hdrs)] * len(hdrs)
    t = ax.table(
        cellText=rows, colLabels=hdrs,
        cellLoc='center', loc='upper center', colWidths=cw,
    )
    t.auto_set_font_size(False)
    t.set_fontsize(9.0)
    t.scale(1, 1.75)
    # header row
    for j in range(len(hdrs)):
        cell = t[0, j]
        cell.set_facecolor(C['th_bg'])
        cell.set_text_props(color=C['bright'], fontweight='bold')
        cell.set_edgecolor(C['hdr'])
        cell.set_linewidth(0)
    # data rows
    for i, row in enumerate(rows):
        bg = C['row_e'] if i % 2 == 0 else C['row_o']
        for j in range(len(hdrs)):
            cell = t[i + 1, j]
            cell.set_facecolor(bg)
            cell.set_edgecolor(C['hdr'])
            cell.set_text_props(color=C['text'])
            cell.set_linewidth(0)
    if title:
        ax.set_title(title, color=C['bright'], fontsize=11,
                     fontweight='bold', pad=4, loc='left')

def textblock(ax, lines, y0=0.97, dy=0.085):
    """
    Render a block of text lines into an axes.
    lines = list of (text, color, fontsize, bold)
    """
    ax.set_facecolor(C['panel'])
    ax.axis('off')
    y = y0
    for text, color, fsize, bold in lines:
        if not text:
            y -= dy * 0.4
            continue
        ax.text(0.01, y, text,
                transform=ax.transAxes,
                color=color, fontsize=fsize,
                fontweight='bold' if bold else 'normal',
                va='top', fontfamily='monospace')
        y -= dy

def embed(fig, parent, toolbar=True):
    """Embed a Figure into a tkinter frame."""
    canvas = FigureCanvasTkAgg(fig, master=parent)
    canvas.draw()
    if toolbar:
        tb = NavigationToolbar2Tk(canvas, parent)
        tb.configure(background=C['hdr'])
        tb.update()
        tb.pack(side=tk.BOTTOM, fill=tk.X)
    canvas.get_tk_widget().configure(bg=C['bg'])
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    return canvas

# ══════════════════════════════════════════════════════════════════
#  DATA LOADING
# ══════════════════════════════════════════════════════════════════
def load_data(path):
    df = pd.read_csv(path)
    # Ensure all key columns are strings before creating the group
    df['CellLine'] = df['CellLine'].astype(str)
    df['Condition'] = df['Condition'].astype(str)

    # Filter out any rows where AspectRatio is missing
    df = df.dropna(subset=['AspectRatio'])

    df['Condition'] = pd.Categorical(
        df['Condition'], categories=COND_ORDER, ordered=True
    )
    df['Sample'] = df['Sample'].astype(str)
    df['group'] = df['CellLine'] + ' | ' + df['Condition'].astype(str)
    return df

# ══════════════════════════════════════════════════════════════════
#  STATISTICS  (all computed once, silently, before the UI opens)
# ══════════════════════════════════════════════════════════════════
def run_stats(df):
    R = {}
    R['cl']   = sorted(df['CellLine'].unique())
    R['cond'] = list(COND_ORDER)

    # ── Two-Way ANOVA ─────────────────────────────────────────────
    import patsy
    model = smf.ols(
        'Q("AspectRatio") ~ Q("CellLine") * Q("Condition")',
        data=df
    ).fit()
    atbl = anova_lm(model, typ=2)
    R['model'] = model
    R['atbl']  = atbl

    k_cl   = 'Q("CellLine")'
    k_cond = 'Q("Condition")'
    k_int  = 'Q("CellLine"):Q("Condition")'   # robust to column order
    k_res  = 'Residual'

    ss = dict(
        cl   = atbl.loc[k_cl,   'sum_sq'],
        cond = atbl.loc[k_cond, 'sum_sq'],
        int  = atbl.loc[k_int,  'sum_sq'],
        err  = atbl.loc[k_res,  'sum_sq'],
    )
    ss['tot'] = sum(ss.values())

    df_d = dict(
        cl   = atbl.loc[k_cl,   'df'],
        cond = atbl.loc[k_cond, 'df'],
        int  = atbl.loc[k_int,  'df'],
        err  = atbl.loc[k_res,  'df'],
    )
    Fv = dict(cl=atbl.loc[k_cl,'F'], cond=atbl.loc[k_cond,'F'], int=atbl.loc[k_int,'F'])
    pv = dict(
        cl   = atbl.loc[k_cl,   'PR(>F)'],
        cond = atbl.loc[k_cond, 'PR(>F)'],
        int  = atbl.loc[k_int,  'PR(>F)'],
    )
    R['ss'] = ss;  R['df'] = df_d;  R['F'] = Fv;  R['p'] = pv

    # Effect sizes
    R['eta2']  = {k: ss[k] / ss['tot']           for k in ['cl','cond','int']}
    R['peta2'] = {k: ss[k] / (ss[k] + ss['err']) for k in ['cl','cond','int']}

    # ── Normality (Lilliefors per group) ──────────────────────────
    norm = []
    for cl in R['cl']:
        for cd in R['cond']:
            d = df.loc[
                (df['CellLine'] == cl) & (df['Condition'] == cd),
                'AspectRatio'
            ].values

            # Add this check to prevent the crash
            if len(d) >= 4:
                stat, p = sm_lilliefors(d)
                norm.append({'label': f'{cl}  |  {cd}', 'n': len(d),
                             'stat': stat, 'p': p, 'data': d, 'cl': cl})
            else:
                # Provide a placeholder for empty/small groups
                print(f"Skipping normality for {cl}|{cd}: only {len(d)} cells found.")
                norm.append({'label': f'{cl}  |  {cd}', 'n': len(d),
                             'stat': np.nan, 'p': np.nan, 'data': d, 'cl': cl})
    R['norm'] = norm

    # ── Levene's test ─────────────────────────────────────────────
    # Get unique groups that actually have data
    active_groups = [g for g in df['group'].unique() if pd.notna(g) and len(df[df['group'] == g]) > 0]
    active_groups.sort()  # This won't crash now because everything is a string

    grps = [
        df.loc[df['group'] == g, 'AspectRatio'].values
        for g in active_groups
    ]

    if len(grps) > 1:
        lstat, lp = stats.levene(*grps, center='median')
    else:
        lstat, lp = np.nan, np.nan

    R['levene'] = dict(stat=lstat, p=lp)
    R['grp_stats'] = [
        {'group': g, 'n': len(d), 'mean': d.mean(), 'std': d.std(), 'var': d.var()}
        for g, d in zip(active_groups, grps)
    ]

    # ── Per-group means & SEMs ────────────────────────────────────
    nc, nd = len(R['cl']), len(R['cond'])
    means = np.zeros((nc, nd));  sems = np.zeros((nc, nd))
    for i, cl in enumerate(R['cl']):
        for j, cd in enumerate(R['cond']):
            d = df.loc[
                (df['CellLine'] == cl) & (df['Condition'] == cd),
                'AspectRatio'
            ].values
            means[i, j] = d.mean()
            sems[i, j]  = d.std() / np.sqrt(len(d))
    R['means'] = means;  R['sems'] = sems

    # ── Tukey-Kramer post-hoc ─────────────────────────────────────
    # Create a mask to only include rows with a valid group name and AspectRatio
    mask = df['group'].notna() & (df['group'].astype(str) != 'nan')
    clean_df = df[mask]

    # Run Tukey only on the cleaned data
    R['tukey'] = pairwise_tukeyhsd(clean_df['AspectRatio'], clean_df['group'].astype(str))

    # ── Linear Mixed-Effects Model ────────────────────────────────
    try:
        lme = smf.mixedlm(
            'AspectRatio ~ C(CellLine) * C(Condition)',
            data=df, groups=df['Sample']
        ).fit(reml=True)
        R['lme'] = lme
    except Exception as e:
        R['lme'] = None;  R['lme_err'] = str(e)

    return R


# ══════════════════════════════════════════════════════════════════
#  APPLICATION
# ══════════════════════════════════════════════════════════════════
class App(tk.Tk):
    def __init__(self, df, R):
        super().__init__()
        self.df = df
        self.R  = R

        self.title('TIBD  Statistical Analysis Suite')
        self.configure(bg=C['hdr'])
        try:    self.state('zoomed')
        except: self.geometry('1420x900')

        self._build_header()
        self._build_notebook()

    # ── chrome ───────────────────────────────────────────────────
    def _build_header(self):
        hdr = tk.Frame(self, bg=C['hdr'], height=62)
        hdr.pack(fill=tk.X)
        hdr.pack_propagate(False)

        tk.Label(
            hdr,
            text='  ◉  TIBD  Statistical Analysis Suite',
            bg=C['hdr'], fg=C['bright'],
            font=('Courier New', 15, 'bold'),
            padx=18
        ).pack(side=tk.LEFT, pady=16)

        info = (
            f"{len(self.df):,} cells  ·  "
            f"{len(self.df['group'].unique())} groups  ·  "
            f"{len(self.R['cl'])} cell lines  ·  "
            f"3 scaffold conditions"
        )
        tk.Label(
            hdr, text=info,
            bg=C['hdr'], fg=C['muted'],
            font=('Courier New', 9),
            padx=18
        ).pack(side=tk.RIGHT, pady=16)

    def _build_notebook(self):
        sty = ttk.Style(self)
        sty.theme_use('clam')
        sty.configure('Dark.TNotebook',
                       background=C['hdr'], borderwidth=0, tabmargins=[0, 0, 0, 0])
        sty.configure('Dark.TNotebook.Tab',
                       background=C['card'], foreground=C['muted'],
                       padding=[16, 9], font=('Courier New', 9, 'bold'),
                       borderwidth=0)
        sty.map('Dark.TNotebook.Tab',
                background=[('selected', C['blue'])],
                foreground=[('selected', C['bright'])])

        nb = ttk.Notebook(self, style='Dark.TNotebook')
        nb.pack(fill=tk.BOTH, expand=True, padx=0, pady=0)

        for label, fn in [
            ('① ANOVA',          self._t_anova),
            ('② Normality',      self._t_normality),
            ("③ Levene's",       self._t_levene),
            ('④ Interaction',    self._t_interaction),
            ('⑤ Post-Hoc',      self._t_posthoc),
            ('⑥ Effect Sizes',   self._t_effects),
            ('⑦ Mixed-Effects',  self._t_lme),
            ('⑧ Box Chart',      self._t_boxchart),
            ('⑨ Violins',        self._t_violin),
            ('⑩ Mean ± SEM',     self._t_meansem),
        ]:
            fr = tk.Frame(nb, bg=C['panel'])
            nb.add(fr, text=f'  {label}  ')
            fn(fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ①  —  TWO-WAY ANOVA
    # ══════════════════════════════════════════════════════════════
    def _t_anova(self, fr):
        R   = self.R
        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])

        ax_t = fig.add_axes([0.03, 0.44, 0.94, 0.53])
        ax_x = fig.add_axes([0.03, 0.02, 0.94, 0.37])

        ms = lambda s, d: s / d
        total_df = int(round(
            R['df']['cl'] + R['df']['cond'] + R['df']['int'] + R['df']['err']
        ))
        rows = [
            ['CellLine',
             f"{R['ss']['cl']:.3f}", f"{int(round(R['df']['cl']))}",
             f"{ms(R['ss']['cl'], R['df']['cl']):.3f}",
             f"{R['F']['cl']:.2f}", f"{R['p']['cl']:.4e}", sig(R['p']['cl'])],
            ['Condition',
             f"{R['ss']['cond']:.3f}", f"{int(round(R['df']['cond']))}",
             f"{ms(R['ss']['cond'], R['df']['cond']):.3f}",
             f"{R['F']['cond']:.2f}", f"{R['p']['cond']:.4e}", sig(R['p']['cond'])],
            ['CellLine × Condition',
             f"{R['ss']['int']:.3f}", f"{int(round(R['df']['int']))}",
             f"{ms(R['ss']['int'], R['df']['int']):.3f}",
             f"{R['F']['int']:.2f}", f"{R['p']['int']:.4e}", sig(R['p']['int'])],
            ['Error',
             f"{R['ss']['err']:.3f}", f"{int(round(R['df']['err']))}",
             f"{ms(R['ss']['err'], R['df']['err']):.3f}", '—', '—', '—'],
            ['Total',
             f"{R['ss']['tot']:.3f}", f"{total_df}", '—', '—', '—', '—'],
        ]
        draw_table(
            ax_t,
            ['Source', 'Sum Sq.', 'd.f.', 'Mean Sq.', 'F', 'Prob > F', 'Sig.'],
            rows,
            cw=[0.24, 0.12, 0.07, 0.12, 0.10, 0.17, 0.08],
            title='Two-Way ANOVA  —  Aspect Ratio ~ CellLine × Condition',
        )

        p_cl, p_cd, p_int = R['p']['cl'], R['p']['cond'], R['p']['int']
        textblock(ax_x, [
            ('━━  RESULTS  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('', '', 0, False),
            (f'  CellLine main effect:              p = {p_cl:.4e}    {sig(p_cl)}', C['blue'],   10.5, False),
            (f'  Condition main effect:             p = {p_cd:.4e}    {sig(p_cd)}', C['blue'],   10.5, False),
            (f'  CellLine × Condition interaction:  p = {p_int:.4e}   {sig(p_int)}   ← significant interaction detected', C['red'], 10.5, False),
            ('', '', 0, False),
            ('━━  INTERPRETATION  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('', '', 0, False),
            ('  ⚠  A significant interaction means the scaffold effect DEPENDS on the cell line.', C['red'],   10, False),
            ('     Main effects CANNOT be interpreted in isolation.', C['text'],  10, False),
            ('     → See Tab ④ for the interaction plot and Tab ⑤ for all pairwise comparisons.', C['muted'], 9.5, False),
            ('', '', 0, False),
            ('  Significance key:  *** p < 0.001   ** p < 0.01   * p < 0.05   ns p ≥ 0.05', C['muted'],  9.5, False),
        ])
        embed(fig, fr, toolbar=False)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ②  —  NORMALITY
    # ══════════════════════════════════════════════════════════════
    def _t_normality(self, fr):
        R    = self.R
        norm = R['norm']
        nCL  = len(R['cl']);  nCd = len(R['cond'])

        fig = plt.Figure(figsize=(14, 9.2), facecolor=C['panel'])
        fig.suptitle(
            'Normality Check  —  Lilliefors Test + QQ Plots',
            color=C['bright'], fontsize=14, fontweight='bold', y=0.995, x=0.02, ha='left'
        )
        fig.text(
            0.02, 0.965,
            '⚠  With N > 10,000 cells, Lilliefors nearly always rejects.  '
            'QQ plots are the true diagnostic — look for heavy tails or bimodality, not statistical perfection.',
            color=C['amber'], fontsize=9.5,
        )

        # QQ grid: nCL rows × nCd cols
        qqW = 0.26; qqH = 0.26
        x_starts = [0.04 + j * (qqW + 0.04) for j in range(nCd)]
        y_starts = [0.64, 0.33]  # row 0 (top), row 1

        idx = 0
        for i in range(nCL):
            for j in range(nCd):
                ax = fig.add_axes([x_starts[j], y_starts[i], qqW, qqH])
                nd  = norm[idx]
                (osm, osr), (slope, intercept, _) = probplot(nd['data'], dist='norm')
                ax.plot(osm, osr, 'o', color=PAL[i], ms=1.3, alpha=0.30, rasterized=True)
                ax.plot(osm, slope * np.array(osm) + intercept,
                        '-', color=C['bright'], lw=1.6, alpha=0.8)
                polish(ax, title=f"{nd['label']}  (n={nd['n']:,})",
                       xl='Theoretical Quantiles', yl='Sample Quantiles', fs=8.5)
                vc = C['red'] if nd['p'] < 0.05 else C['green']
                vt = 'Non-Normal ✗' if nd['p'] < 0.05 else 'Normal ✓'
                ax.text(0.97, 0.04, vt, transform=ax.transAxes,
                        ha='right', va='bottom', color=vc, fontsize=8, fontweight='bold')
                idx += 1

        # Lilliefors summary table
        ax_tbl = fig.add_axes([0.03, 0.02, 0.94, 0.25])
        tbl_rows = [
            [nd['label'], f"{nd['n']:,}", f"{nd['stat']:.5f}",
             f"{nd['p']:.3e}", 'Non-Normal ✗' if nd['p'] < 0.05 else 'Normal ✓']
            for nd in norm
        ]
        draw_table(
            ax_tbl,
            ['Group', 'n', 'Test Statistic', 'p-value', 'Result'],
            tbl_rows,
            cw=[0.37, 0.10, 0.18, 0.18, 0.17],
            title='Lilliefors Test Summary',
        )
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ③  —  LEVENE'S TEST
    # ══════════════════════════════════════════════════════════════
    def _t_levene(self, fr):
        R   = self.R
        lev = R['levene']
        gs  = R['grp_stats']

        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])
        fig.suptitle(
            "Levene's Test  —  Homogeneity of Variance",
            color=C['bright'], fontsize=14, fontweight='bold', y=0.98, x=0.02, ha='left'
        )

        # Variance bar chart
        ax_v = fig.add_axes([0.04, 0.38, 0.50, 0.52])
        labels    = [g['group'] for g in gs]
        variances = [g['var']   for g in gs]
        bar_cols  = [PAL[0] if R['cl'][0] in g else PAL[1] for g in labels]
        bars = ax_v.barh(range(len(labels)), variances,
                         color=bar_cols, alpha=0.82,
                         edgecolor=C['panel'], linewidth=2)
        ax_v.set_yticks(range(len(labels)))
        ax_v.set_yticklabels(labels, fontsize=9.5, color=C['text'])
        ax_v.invert_yaxis()
        for bar, val in zip(bars, variances):
            ax_v.text(val + max(variances)*0.01, bar.get_y() + bar.get_height()/2,
                      f'{val:.5f}', va='center', fontsize=8.5, color=C['muted'])
        polish(ax_v, title='Variance per Group', xl='Variance')

        # Descriptive stats table
        ax_t = fig.add_axes([0.58, 0.38, 0.40, 0.52])
        tbl_rows = [
            [g['group'], f"{g['n']:,}", f"{g['mean']:.4f}",
             f"{g['std']:.4f}", f"{g['var']:.5f}"]
            for g in gs
        ]
        draw_table(
            ax_t,
            ['Group', 'n', 'Mean', 'Std Dev', 'Variance'],
            tbl_rows,
            cw=[0.36, 0.12, 0.16, 0.16, 0.20],
            title='Group Descriptive Statistics',
        )

        # Result + rationale
        ax_x = fig.add_axes([0.04, 0.02, 0.93, 0.31])
        vc = C['red'] if lev['p'] < 0.05 else C['green']
        vt = 'HETEROGENEOUS  (p < 0.05)' if lev['p'] < 0.05 else 'HOMOGENEOUS  (p ≥ 0.05)'
        textblock(ax_x, [
            (f"  F = {lev['stat']:.4f}   ·   p = {lev['p']:.4e}   ·   {sig(lev['p'])}", C['text'],  11, False),
            (f"  Verdict:  {vt}", vc, 12, True),
            ('', '', 0, False),
            ("━━  WHY LEVENE'S AND NOT BARTLETT'S?  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━", C['bright'], 10, True),
            ("  Levene's uses absolute deviations from group medians — robust to non-normality.", C['text'], 9.5, False),
            ("  Bartlett's assumes normality (violated here) → would give unreliable results.", C['text'], 9.5, False),
            ('', '', 0, False),
            ('━━  ANOVA IMPLICATION  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('  Heterogeneous variances reduce ANOVA F-test reliability.', C['text'], 9.5, False),
            ('  The Linear Mixed-Effects Model (Tab ⑦) provides more valid inference in this case.', C['text'], 9.5, False),
        ], dy=0.11)
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ④  —  INTERACTION PLOT
    # ══════════════════════════════════════════════════════════════
    def _t_interaction(self, fr):
        R = self.R
        x = np.arange(len(R['cond']))

        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])
        ax  = fig.add_axes([0.07, 0.11, 0.60, 0.82])

        for i, cl in enumerate(R['cl']):
            m, s = R['means'][i], R['sems'][i]
            ax.fill_between(x, m - s, m + s, color=PAL[i], alpha=0.15)
            ax.errorbar(x, m, yerr=s,
                        fmt='-o', color=PAL[i], label=cl,
                        lw=2.8, ms=10, capsize=7, capthick=2.2,
                        markerfacecolor=PAL[i], markeredgecolor=C['panel'],
                        markeredgewidth=2, elinewidth=2, zorder=5)

        ax.set_xticks(x)
        ax.set_xticklabels(R['cond'], fontsize=12)
        polish(ax,
               title='Interaction:  CellLine × Scaffold Condition',
               xl='Scaffold Condition',
               yl='Mean Aspect Ratio ± SEM')
        ax.legend(title='Cell Line', fontsize=11, title_fontsize=10.5)

        # Info sidebar
        ax_i = fig.add_axes([0.71, 0.11, 0.26, 0.82])
        ax_i.set_facecolor(C['card'])
        ax_i.axis('off')
        p_int = R['p']['int']
        textblock(ax_i, [
            ('RESULT', C['bright'], 11, True),
            (f"p = {p_int:.4e}", C['red'], 11, False),
            (f"  {sig(p_int)}", C['red'], 18, True),
            ('', '', 0, False),
            ('INTERPRETATION', C['bright'], 10, True),
            ('', '', 0, False),
            ('Non-parallel lines confirm', C['text'], 9.5, False),
            ('a significant interaction.', C['text'], 9.5, False),
            ('', '', 0, False),
            ('The scaffold effect', C['text'], 9.5, False),
            ('DEPENDS on cell line.', C['red'],  9.5, True),
            ('', '', 0, False),
            ('Main effects cannot be', C['text'], 9.5, False),
            ('interpreted in isolation.', C['text'], 9.5, False),
            ('', '', 0, False),
            ('→ Tab ⑤ for pairwise', C['muted'], 9, False),
            ('  significance testing', C['muted'], 9, False),
        ], dy=0.062)
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ⑤  —  TUKEY-KRAMER POST-HOC
    # ══════════════════════════════════════════════════════════════
    def _t_posthoc(self, fr):
        R     = self.R
        tukey = R['tukey']
        raw   = tukey.summary().data[1:]   # skip header row

        g1   = [str(r[0]) for r in raw]
        g2   = [str(r[1]) for r in raw]
        md   = [float(r[2]) for r in raw]
        padj = [float(r[3]) for r in raw]
        lo   = [float(r[4]) for r in raw]
        hi   = [float(r[5]) for r in raw]
        rej  = [str(r[6]).lower() == 'true' for r in raw]
        n    = len(g1)

        fig = plt.Figure(figsize=(14, 10.5), facecolor=C['panel'])
        fig.suptitle(
            'Tukey-Kramer Post-Hoc  —  All Pairwise Comparisons',
            color=C['bright'], fontsize=14, fontweight='bold', y=0.995, x=0.02, ha='left'
        )

        # Forest plot
        ax_f  = fig.add_axes([0.03, 0.05, 0.55, 0.92])
        y_pos = np.arange(n, 0, -1)

        for i in range(n):
            col   = C['red']  if rej[i] else C['muted']
            alpha = 1.0 if rej[i] else 0.45
            ax_f.plot([lo[i], hi[i]], [y_pos[i]] * 2,
                      '-', color=col, lw=2.6, alpha=alpha, zorder=4)
            ax_f.plot(md[i], y_pos[i], 'o',
                      color=col, ms=7, alpha=alpha,
                      markerfacecolor=col, markeredgecolor=C['panel'],
                      mew=1.5, zorder=5)

        ax_f.axvline(0, color=C['bright'], lw=1.4, ls='--', alpha=0.5, zorder=3)
        lbl = [
            f"{g1[i].split('|')[-1].strip()} vs {g2[i].split('|')[-1].strip()}"
            for i in range(n)
        ]
        ax_f.set_yticks(y_pos)
        ax_f.set_yticklabels(lbl, fontsize=8.5)
        polish(ax_f,
               title='Comparison Intervals  (Red = Significant)',
               xl='Mean Difference  (95% Tukey-Kramer CI)')
        ax_f.legend(handles=[
            mpatches.Patch(color=C['red'],  label='Significant  p < 0.05'),
            mpatches.Patch(color=C['muted'],label='Not significant'),
        ], loc='lower right', fontsize=9)

        # Results table
        ax_t = fig.add_axes([0.61, 0.05, 0.38, 0.92])
        tbl  = [
            [g1[i].split('|')[-1].strip(),
             g2[i].split('|')[-1].strip(),
             f'{md[i]:+.4f}', f'{padj[i]:.3e}', sig(padj[i])]
            for i in range(n)
        ]
        draw_table(
            ax_t,
            ['Group 1', 'Group 2', 'Δ Mean', 'Adj. p', 'Sig.'],
            tbl,
            cw=[0.25, 0.25, 0.18, 0.22, 0.10],
            title='All Pairwise Results',
        )
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ⑥  —  EFFECT SIZES
    # ══════════════════════════════════════════════════════════════
    def _t_effects(self, fr):
        R       = self.R
        sources = ['CellLine', 'Condition', 'Interaction']
        ev      = [R['eta2']['cl'],  R['eta2']['cond'],  R['eta2']['int']]
        pev     = [R['peta2']['cl'], R['peta2']['cond'], R['peta2']['int']]

        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])
        fig.suptitle(
            'Effect Sizes  —  η²  (Eta²)  and  Partial η²',
            color=C['bright'], fontsize=14, fontweight='bold', y=0.98, x=0.02, ha='left'
        )

        ax  = fig.add_axes([0.07, 0.30, 0.55, 0.60])
        x   = np.arange(3);  w = 0.30
        b1  = ax.bar(x - w/2, ev,  w, label='η²',
                     color=C['blue'],   alpha=0.85, edgecolor=C['panel'], zorder=4)
        b2  = ax.bar(x + w/2, pev, w, label='Partial η²',
                     color=C['orange'], alpha=0.85, edgecolor=C['panel'], zorder=4)
        for bar, val in zip(list(b1) + list(b2), ev + pev):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.0004,
                    f'{val:.4f}', ha='center', va='bottom', fontsize=9, color=C['text'])
        for thr, lbl_t, ls in [
            (0.01, 'Small',   ':'),
            (0.06, 'Medium', '--'),
            (0.14, 'Large',   '-'),
        ]:
            ax.axhline(thr, color=C['muted'], lw=1.0, ls=ls, alpha=0.6)
            ax.text(2.66, thr + 0.0006, lbl_t, fontsize=8, color=C['muted'], va='bottom')
        ax.set_xticks(x)
        ax.set_xticklabels(sources, fontsize=11)
        polish(ax, title='Effect Size by Source', yl='Effect Size')
        ax.legend(fontsize=11, framealpha=0.85)

        # Summary table
        ax_t = fig.add_axes([0.66, 0.30, 0.32, 0.60])
        draw_table(
            ax_t,
            ['Source', 'η²', 'Partial η²', 'Magnitude'],
            [[s, f'{e:.4f}', f'{pe:.4f}', emag(pe)]
             for s, e, pe in zip(sources, ev, pev)],
            cw=[0.30, 0.17, 0.22, 0.31],
            title='Summary',
        )

        # Reference text
        ax_x = fig.add_axes([0.04, 0.02, 0.93, 0.23])
        textblock(ax_x, [
            ('━━  REFERENCE THRESHOLDS (Cohen 1988)  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('  Negligible: pη² < 0.01   ·   Small: pη² ≥ 0.01   ·   Medium: pη² ≥ 0.06   ·   Large: pη² ≥ 0.14', C['text'], 9.5, False),
            ('', '', 0, False),
            ('━━  FORMULAS  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('  η²  = SS_effect / SS_total                  →  raw share of total variance', C['text'], 9.5, False),
            ('  pη² = SS_effect / (SS_effect + SS_error)    →  preferred for factorial designs; isolates each factor', C['text'], 9.5, False),
        ], dy=0.18)
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ⑦  —  LINEAR MIXED-EFFECTS MODEL
    # ══════════════════════════════════════════════════════════════
    def _t_lme(self, fr):
        R   = self.R
        lme = R.get('lme')
        fig = plt.Figure(figsize=(14, 9), facecolor=C['panel'])
        fig.suptitle(
            'Linear Mixed-Effects  —  AspectRatio ~ CellLine × Condition + (1 | Sample)',
            color=C['bright'], fontsize=13, fontweight='bold', y=0.995, x=0.02, ha='left'
        )

        if lme is None:
            ax = fig.add_axes([0.1, 0.3, 0.8, 0.4])
            ax.axis('off')
            ax.text(0.5, 0.5, f"Model fitting failed:\n{R.get('lme_err', 'unknown')}",
                    ha='center', va='center', color=C['red'], fontsize=12)
            embed(fig, fr, toolbar=False)
            return

        # Fixed effects table
        fe_rows = []
        for pn in lme.params.index:
            pval = lme.pvalues[pn]
            fe_rows.append([
                pn[:52], f'{lme.params[pn]:.4f}', f'{lme.bse[pn]:.4f}',
                f'{lme.tvalues[pn]:.3f}', f'{pval:.4e}', sig(pval),
            ])
        ax_fe = fig.add_axes([0.03, 0.60, 0.94, 0.33])
        draw_table(
            ax_fe,
            ['Term', 'Estimate', 'SE', 't-Stat', 'p-Value', 'Sig.'],
            fe_rows,
            cw=[0.44, 0.12, 0.10, 0.10, 0.16, 0.08],
            title='Fixed-Effects Coefficients',
        )

        # Random effects
        try:
            re_var = float(lme.cov_re.iloc[0, 0])
        except Exception:
            try:
                re_var = float(lme.vcomp[0])
            except Exception:
                re_var = float('nan')

        ax_re = fig.add_axes([0.03, 0.45, 0.94, 0.12])
        draw_table(
            ax_re,
            ['Component', 'Variance', 'Description'],
            [
                ['Sample  (random intercept)', f'{re_var:.6f}',  'Image-level clustering variance'],
                ['Residual',                   f'{lme.scale:.6f}', 'Within-image cell-to-cell variance'],
            ],
            cw=[0.28, 0.18, 0.54],
            title='Random-Effects Variance Components',
        )

        # Rationale
        ax_x = fig.add_axes([0.03, 0.02, 0.94, 0.38])
        textblock(ax_x, [
            ('━━  WHY A MIXED-EFFECTS MODEL?  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('  Each image contributes many cells that share the same imaging conditions and local microenvironment.', C['text'], 9.5, False),
            ('  Those cells are CORRELATED — treating them as independent inflates effective N and produces', C['text'], 9.5, False),
            ('  artificially small p-values (pseudoreplication).', C['red'],  9.5, False),
            ('', '', 0, False),
            ('  The (1|Sample) random intercept absorbs image-level clustering, giving correct degrees of freedom.', C['text'], 9.5, False),
            ('', '', 0, False),
            ('━━  VALIDATION RULE  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━', C['bright'], 10, True),
            ('  ANOVA and LME agree   →  pseudoreplication does not distort conclusions.    Trust either.', C['green'], 9.5, False),
            ('  They disagree         →  pseudoreplication is real.    Trust the LME.', C['red'],   9.5, False),
        ], dy=0.115)
        embed(fig, fr, toolbar=False)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ⑧  —  BOX CHART
    # ══════════════════════════════════════════════════════════════
    def _t_boxchart(self, fr):
        R = self.R;  df = self.df
        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])
        ax  = fig.add_axes([0.07, 0.11, 0.77, 0.83])

        medianprops = dict(color=C['bright'], linewidth=2.5)
        whiskerprops = dict(color=C['muted'], linewidth=1.5, linestyle='--')
        capprops = dict(color=C['muted'], linewidth=1.5)
        boxprops = dict(linewidth=1.5)

        sns.boxplot(
            data=df, x='Condition', y='AspectRatio', hue='CellLine',
            order=R['cond'], palette=PAL,
            width=0.60, linewidth=1.3, fliersize=0, ax=ax,
            medianprops=medianprops, whiskerprops=whiskerprops,
            capprops=capprops, boxprops=boxprops,
        )
        polish(
            ax,
            title='Cell Elongation Distribution  —  Aspect Ratio by Condition & Cell Line',
            xl='Scaffold Condition',
            yl='Aspect Ratio  (Length / Width)',
        )
        ax.legend(title='Cell Line', fontsize=11, title_fontsize=10.5, framealpha=0.85)
        ax.set_ylim(0.85, 4.3)

        # Sample size annotations
        for i, cd in enumerate(R['cond']):
            for j, cl in enumerate(R['cl']):
                n = len(df[(df['Condition'] == cd) & (df['CellLine'] == cl)])
                ax.text(i + [-0.20, 0.20][j], 0.92,
                        f'n={n:,}', ha='center', va='bottom',
                        fontsize=7.5, color=PAL[j], fontweight='bold')
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ⑨  —  VIOLIN PLOTS
    # ══════════════════════════════════════════════════════════════
    def _t_violin(self, fr):
        R = self.R;  df = self.df
        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])
        ax  = fig.add_axes([0.07, 0.11, 0.77, 0.83])

        vio_kw = dict(
            data=df, x='Condition', y='AspectRatio', hue='CellLine',
            order=R['cond'], palette=PAL,
            inner='box', cut=0, linewidth=1.1, ax=ax,
        )
        try:
            sns.violinplot(**vio_kw, density_norm='width')   # seaborn ≥ 0.12
        except TypeError:
            sns.violinplot(**vio_kw, scale='width')           # seaborn < 0.12

        polish(
            ax,
            title='Aspect Ratio Distributions  —  KDE Violin Plot',
            xl='Scaffold Condition',
            yl='Aspect Ratio',
        )
        ax.legend(title='Cell Line', fontsize=11, title_fontsize=10.5, framealpha=0.85)
        ax.set_ylim(0.8, 5.0)
        embed(fig, fr)

    # ══════════════════════════════════════════════════════════════
    #  TAB  ⑩  —  MEAN ± SEM BAR CHART
    # ══════════════════════════════════════════════════════════════
    def _t_meansem(self, fr):
        R = self.R
        means, sems = R['means'], R['sems']

        fig = plt.Figure(figsize=(14, 8.5), facecolor=C['panel'])
        ax  = fig.add_axes([0.07, 0.11, 0.77, 0.83])
        x   = np.arange(len(R['cond']));  w = 0.30
        off = [-w / 2, w / 2]

        for i, cl in enumerate(R['cl']):
            ax.bar(x + off[i], means[i], w,
                   label=cl, color=PAL[i], alpha=0.85,
                   edgecolor=C['panel'], linewidth=2, zorder=4)
            ax.errorbar(x + off[i], means[i], yerr=sems[i],
                        fmt='none', color=C['bright'],
                        lw=2.2, capsize=8, capthick=2.2, zorder=5)
            for xi, (m, s) in enumerate(zip(means[i], sems[i])):
                ax.text(xi + off[i], m + s + 0.022, f'{m:.3f}',
                        ha='center', va='bottom',
                        fontsize=8.5, color=PAL[i], fontweight='bold')

        ax.set_xticks(x)
        ax.set_xticklabels(R['cond'], fontsize=12)
        polish(
            ax,
            title='Group Means with Standard Error of the Mean',
            xl='Scaffold Condition',
            yl='Mean Aspect Ratio ± SEM',
        )
        ax.legend(title='Cell Line', fontsize=11, title_fontsize=10.5, framealpha=0.85)
        # Find the maximum mean while ignoring the NaN values from empty groups
        y_max = np.nanmax(R['means']) if not np.all(np.isnan(R['means'])) else 1.0
        ax.set_ylim(0, y_max * 1.32)
        embed(fig, fr)


# ══════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ══════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    # allow command-line path override
    path = sys.argv[1] if len(sys.argv) > 1 else CSV_PATH

    # fall back to file-picker if path doesn't exist
    if not os.path.exists(path):
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename(
            title='Select BioRep CSV file',
            filetypes=[('CSV files', '*.csv'), ('All files', '*.*')],
        )
        root.destroy()
        if not path:
            sys.exit(0)

    print('Loading data ...')
    df = load_data(path)
    print(f'  {len(df):,} cells loaded across {len(df["group"].unique())} groups.')
    print('Running statistics (this may take a moment for the mixed model) ...')
    R  = run_stats(df)
    print('  Done.  Launching window ...')

    App(df, R).mainloop()