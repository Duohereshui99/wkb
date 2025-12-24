#!/usr/bin/env python3
"""
Generate PRC-format figures - Fixed version with no overlap
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

COLUMN_WIDTH = 3.4

plt.rcParams.update({
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
})

# Data
actinides = {
    'Ra226': {'Z': 88, 'N': 138, 'P_O0': 0.219, 'P_O2': 0.125},
    'Th230': {'Z': 90, 'N': 140, 'P_O0': 0.282, 'P_O2': 0.159},
    'Th232': {'Z': 90, 'N': 142, 'P_O0': 0.419, 'P_O2': 0.235},
    'U234':  {'Z': 92, 'N': 142, 'P_O0': 0.257, 'P_O2': 0.145},
    'U236':  {'Z': 92, 'N': 144, 'P_O0': 0.312, 'P_O2': 0.175},
    'U238':  {'Z': 92, 'N': 146, 'P_O0': 0.422, 'P_O2': 0.237},
    'Pu238': {'Z': 94, 'N': 144, 'P_O0': 0.159, 'P_O2': 0.091},
    'Pu240': {'Z': 94, 'N': 146, 'P_O0': 0.215, 'P_O2': 0.122},
    'Pu242': {'Z': 94, 'N': 148, 'P_O0': 0.218, 'P_O2': 0.123},
    'Pu244': {'Z': 94, 'N': 150, 'P_O0': 0.190, 'P_O2': 0.107},
    'Am241': {'Z': 95, 'N': 146, 'P_O0': 0.052, 'P_O2': 0.030},
    'Am243': {'Z': 95, 'N': 148, 'P_O0': 0.042, 'P_O2': 0.024},
    'Cm244': {'Z': 96, 'N': 148, 'P_O0': 0.119, 'P_O2': 0.067},
    'Cm246': {'Z': 96, 'N': 150, 'P_O0': 0.120, 'P_O2': 0.068},
    'Cm248': {'Z': 96, 'N': 152, 'P_O0': 0.150, 'P_O2': 0.084},
}

superheavy = {
    'Ds275': {'Z': 110, 'N': 165, 'P_O0': 0.0032, 'P_O2': 0.0029},
    'Hs271': {'Z': 108, 'N': 163, 'P_O0': 0.0260, 'P_O2': 0.0235},
    'Sg267': {'Z': 106, 'N': 161, 'P_O0': 0.0225, 'P_O2': 0.0202},
    'Ds273': {'Z': 110, 'N': 163, 'P_O0': 0.0137, 'P_O2': 0.0125},
    'Hs269': {'Z': 108, 'N': 161, 'P_O0': 0.0547, 'P_O2': 0.0494},
    'Sg265': {'Z': 106, 'N': 159, 'P_O0': 0.0220, 'P_O2': 0.0199},
    'Rf261': {'Z': 104, 'N': 157, 'P_O0': 0.0366, 'P_O2': 0.0329},
    'Lr251': {'Z': 103, 'N': 148, 'P_O0': 0.0438, 'P_O2': 0.0396},
    'Lr253': {'Z': 103, 'N': 150, 'P_O0': 0.0354, 'P_O2': 0.0319},
    'Lr255': {'Z': 103, 'N': 152, 'P_O0': 0.0934, 'P_O2': 0.0840},
    'Hs272': {'Z': 108, 'N': 164, 'P_O0': 0.0198, 'P_O2': 0.0179},
    'Ds276': {'Z': 110, 'N': 166, 'P_O0': 0.0427, 'P_O2': 0.0390},
    'Cn284': {'Z': 112, 'N': 172, 'P_O0': 0.0755, 'P_O2': 0.0682},
    'Fl286': {'Z': 114, 'N': 172, 'P_O0': 0.0301, 'P_O2': 0.0273},
    'Fl288': {'Z': 114, 'N': 174, 'P_O0': 0.0547, 'P_O2': 0.0495},
    'Lv292': {'Z': 116, 'N': 176, 'P_O0': 0.1213, 'P_O2': 0.1103},
}

# ============================================
# Figure 1: Preformation comparison
# ============================================
def plot_preformation_comparison():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLUMN_WIDTH, 5.2))
    plt.subplots_adjust(hspace=0.32, left=0.17, right=0.96, top=0.94, bottom=0.09)
    
    act_N = [d['N'] for d in actinides.values()]
    act_P0 = [d['P_O0'] for d in actinides.values()]
    act_P2 = [d['P_O2'] for d in actinides.values()]
    
    sh_N = [d['N'] for d in superheavy.values()]
    sh_P0 = [d['P_O0'] for d in superheavy.values()]
    sh_P2 = [d['P_O2'] for d in superheavy.values()]
    
    # Panel (a)
    ax1.scatter(act_N, act_P0, s=35, marker='s', facecolors='none', edgecolors='blue', 
                linewidths=1.2, label='Act. O0', zorder=3)
    ax1.scatter(act_N, act_P2, s=35, marker='o', facecolors='none', edgecolors='red', 
                linewidths=1.2, label='Act. O2', zorder=3)
    ax1.scatter(sh_N, sh_P0, s=35, marker='s', color='blue', alpha=0.6, 
                label='SH O0', zorder=3)
    ax1.scatter(sh_N, sh_P2, s=35, marker='o', color='red', alpha=0.6, 
                label='SH O2', zorder=3)
    
    # Shell closure line for N=162
    ax1.axvline(x=162, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
    ax1.text(163, 0.2, '$N$=162', fontsize=7, ha='left', va='bottom')
    
    ax1.set_xlabel('Neutron number $N$')
    ax1.set_ylabel('$P_\\alpha$')
    ax1.set_ylim(-0.02, 0.48)
    ax1.set_xlim(135, 180)
    # Legend in upper right with frame, moved down to avoid (a) label
    ax1.legend(loc='upper right', frameon=True, ncol=2, fontsize=7,
               columnspacing=0.5, handletextpad=0.2, fancybox=False, edgecolor='black',
               bbox_to_anchor=(0.99, 0.88))
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.text(0.95, 0.92, '(a)', transform=ax1.transAxes, fontsize=10, fontweight='bold', ha='right')
    
    # Panel (b)
    ax2.scatter(act_P0, act_P2, s=35, marker='s', facecolors='none', edgecolors='blue', 
                linewidths=1.2, label='Actinides', zorder=3)
    ax2.scatter(sh_P0, sh_P2, s=35, marker='o', color='red', alpha=0.6, 
                label='Superheavy', zorder=3)
    
    x_line = np.linspace(0, 0.45, 100)
    ax2.plot(x_line, 0.56 * x_line, 'b--', linewidth=1, alpha=0.8, label='0.56')
    ax2.plot(x_line, 0.90 * x_line, 'r:', linewidth=1.2, alpha=0.8, label='0.90')
    ax2.plot([0, 0.45], [0, 0.45], 'k-', linewidth=0.5, alpha=0.3)
    
    ax2.set_xlabel('$P_\\alpha^{\\mathrm{O0}}$')
    ax2.set_ylabel('$P_\\alpha^{\\mathrm{O2}}$')
    ax2.set_xlim(-0.01, 0.46)
    ax2.set_ylim(-0.01, 0.28)
    ax2.legend(loc='lower right', frameon=False, fontsize=7, ncol=2,
               columnspacing=0.5, handletextpad=0.2)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.text(0.95, 0.92, '(b)', transform=ax2.transAxes, fontsize=10, fontweight='bold', ha='right')
    
    plt.savefig('preformation_comparison.pdf', bbox_inches='tight')
    plt.savefig('preformation_comparison.png', bbox_inches='tight')
    print('Saved preformation_comparison.pdf/png')
    plt.close()

# ============================================
# Figure 2: Shell structure
# ============================================
def plot_shell_structure():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COLUMN_WIDTH, 5.2))
    plt.subplots_adjust(hspace=0.32, left=0.17, right=0.96, top=0.94, bottom=0.09)
    
    # Panel (a): Actinides
    act_dN = [d['N'] - 126 for d in actinides.values()]
    act_P0 = [d['P_O0'] for d in actinides.values()]
    act_P2 = [d['P_O2'] for d in actinides.values()]
    
    ax1.scatter(act_dN, act_P0, s=35, marker='s', facecolors='none', edgecolors='blue', 
                linewidths=1.2, label='O0', zorder=3)
    ax1.scatter(act_dN, act_P2, s=35, marker='o', facecolors='none', edgecolors='red', 
                linewidths=1.2, label='O2', zorder=3)
    
    z0 = np.polyfit(act_dN, act_P0, 1)
    z2 = np.polyfit(act_dN, act_P2, 1)
    x_fit = np.linspace(10, 28, 50)
    ax1.plot(x_fit, np.polyval(z0, x_fit), 'b--', linewidth=1, alpha=0.7)
    ax1.plot(x_fit, np.polyval(z2, x_fit), 'r--', linewidth=1, alpha=0.7)
    
    r0 = np.corrcoef(act_dN, act_P0)[0,1]
    r2 = np.corrcoef(act_dN, act_P2)[0,1]
    
    ax1.set_xlabel('$N - 126$')
    ax1.set_ylabel('$P_\\alpha$')
    ax1.set_title('Actinides', fontsize=9, pad=3)
    ax1.legend(loc='upper right', frameon=True, fontsize=8,
               fancybox=False, edgecolor='black')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.text(0.03, 0.92, '(a)', transform=ax1.transAxes, fontsize=10, fontweight='bold')
    
    # Panel (b): Superheavy
    sh_dN = [d['N'] - 162 for d in superheavy.values()]
    sh_P0 = [d['P_O0'] for d in superheavy.values()]
    sh_P2 = [d['P_O2'] for d in superheavy.values()]
    
    ax2.scatter(sh_dN, sh_P0, s=35, marker='s', facecolors='none', edgecolors='blue', 
                linewidths=1.2, label='O0', zorder=3)
    ax2.scatter(sh_dN, sh_P2, s=35, marker='o', facecolors='none', edgecolors='red', 
                linewidths=1.2, label='O2', zorder=3)
    
    z0_sh = np.polyfit(sh_dN, sh_P0, 1)
    z2_sh = np.polyfit(sh_dN, sh_P2, 1)
    x_fit_sh = np.linspace(-8, 16, 50)
    ax2.plot(x_fit_sh, np.polyval(z0_sh, x_fit_sh), 'b--', linewidth=1, alpha=0.7)
    ax2.plot(x_fit_sh, np.polyval(z2_sh, x_fit_sh), 'r--', linewidth=1, alpha=0.7)
    
    r0_sh = np.corrcoef(sh_dN, sh_P0)[0,1]
    r2_sh = np.corrcoef(sh_dN, sh_P2)[0,1]
    
    ax2.set_xlabel('$N - 162$')
    ax2.set_ylabel('$P_\\alpha$')
    ax2.set_title('Superheavy', fontsize=9, pad=3)
    # No legend for (b), same as (a); correlation values in text
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.text(0.03, 0.92, '(b)', transform=ax2.transAxes, fontsize=10, fontweight='bold')
    
    plt.savefig('shell_structure.pdf', bbox_inches='tight')
    plt.savefig('shell_structure.png', bbox_inches='tight')
    print('Saved shell_structure.pdf/png')
    plt.close()

# ============================================
# Figure 3: S2 contribution
# ============================================
def plot_s2_contribution():
    hbar = 197.327
    m_N = 931.494
    e2 = 1.44

    Z_d, A_d = 90, 234
    Z_alpha = 2
    Q = 4.270

    mu = 4 * A_d / (4 + A_d) * m_N
    r0 = 1.17
    R0 = r0 * (A_d**(1/3) + 4**(1/3))
    a = 0.70
    V0 = 50 + 33 * (A_d - 2*Z_d) / (A_d + 4)

    def V_eff(r):
        V_N = -V0 / (1 + np.exp((r - R0) / a))
        V_C = np.where(r < R0,
                       Z_alpha * Z_d * e2 / (2*R0) * (3 - r**2/R0**2),
                       Z_alpha * Z_d * e2 / r)
        return V_N + V_C

    from scipy.optimize import brentq
    r_in = brentq(lambda r: V_eff(r) - Q, 7, 10)
    r_out = brentq(lambda r: V_eff(r) - Q, 50, 80)

    r = np.linspace(r_in + 0.5, r_out - 0.5, 500)
    V = V_eff(r)
    Q_val = 2 * mu / hbar**2 * (Q - V)

    dr = r[1] - r[0]
    Q_prime = np.gradient(Q_val, dr)
    Q_double_prime = np.gradient(Q_prime, dr)

    S0 = np.sqrt(np.abs(Q_val))
    S2 = Q_double_prime / (8 * np.abs(Q_val)**1.5) - 5 * Q_prime**2 / (32 * np.abs(Q_val)**2.5)

    # Larger figure with bigger fonts for PRC double-column readability
    fig, ax = plt.subplots(figsize=(COLUMN_WIDTH, 3.0))
    plt.subplots_adjust(left=0.20, right=0.96, top=0.86, bottom=0.22)

    ax.plot(r, S2/S0 * 100, 'b-', linewidth=1.8)
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    ax.fill_between(r, S2/S0 * 100, 0, where=(S2 < 0), alpha=0.2, color='green')

    ax.set_xlabel('$r$ (fm)', fontsize=14)
    ax.set_ylabel('$S_2/S_0$ (%)', fontsize=14)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlim(r_in - 1, r_out + 1)
    ax.set_ylim(-3.5, 1.5)

    ax.axvline(x=r_in, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axvline(x=r_out, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.text(r_in + 2, 1.0, '$r_{\\mathrm{in}}$', fontsize=12)
    ax.text(r_out - 8, 1.0, '$r_{\\mathrm{out}}$', fontsize=12)

    ax.set_title('$^{238}$U', fontsize=14, pad=6)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    plt.savefig('s2_contribution_U238.pdf', bbox_inches='tight')
    plt.savefig('s2_contribution_U238.png', bbox_inches='tight')
    print('Saved s2_contribution_U238.pdf/png')
    plt.close()

if __name__ == '__main__':
    plot_preformation_comparison()
    plot_shell_structure()
    plot_s2_contribution()
    print('\nAll figures done!')
