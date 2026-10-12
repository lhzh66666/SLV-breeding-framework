# -*- coding: utf-8 -*-
"""
大豆SLV框架时间盲测：完整最终版（含数据清洗、BLUP、ΔS、验证、LOYO、敏感性分析）
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============ 1. 读取数据 ============
print("正在读取数据...")
pheno = pd.read_excel("Table1-4.xlsx", sheet_name="Sheet 1  Phenotype（表型数据）", skiprows=1)
pheno.columns = ["Line_ID", "Year", "Location", "Yield", "Seed_Weight", "Isoflavone"]
pheno["Line_ID"] = pheno["Line_ID"].astype(str)
pheno = pheno.dropna(subset=["Yield", "Isoflavone", "Seed_Weight"])
pheno["Year"] = pheno["Year"].astype(int)

print(f"家系数：{pheno['Line_ID'].nunique()}")
print(f"年份：{sorted(pheno['Year'].unique())}")
print(f"地点：{sorted(pheno['Location'].unique())}")
print(f"总观测数：{len(pheno)}")

# ============ 2. 划分训练集和验证集 ============
train_years = [2006, 2011, 2012, 2014]
valid_year = 2015
train = pheno[pheno["Year"].isin(train_years)].copy()
valid = pheno[pheno["Year"] == valid_year].copy()
print(f"\n训练集观测数：{len(train)}")
print(f"验证集观测数：{len(valid)}")

# ============ 3. 计算BLUP（基于品系均值） ============
def calc_blup(df, traits):
    blup_dict = {}
    for tr in traits:
        df_tr = df[["Line_ID", "Location", tr]].dropna()
        df_tr["centered"] = df_tr[tr] - df_tr.groupby("Location")[tr].transform("mean")
        line_effect = df_tr.groupby("Line_ID")["centered"].mean()
        blup_dict[tr] = line_effect
    return pd.DataFrame(blup_dict)

traits = ["Yield", "Isoflavone", "Seed_Weight"]
blup_train = calc_blup(train, traits)
blup_valid = calc_blup(valid, traits)

common_lines = blup_train.index.intersection(blup_valid.index)
blup_train = blup_train.loc[common_lines]
blup_valid = blup_valid.loc[common_lines]
print(f"\n共有家系数：{len(common_lines)}")

# ============ 4. 计算遗传相关矩阵 ============
r_matrix = blup_train.corr()
print("\n遗传相关系数矩阵：")
print(r_matrix.round(3))

# ============ 5. 计算ΔS ============
omega = np.array([0.20, 0.65, 0.15])
lambda_ = 0.3

delta_T_train = blup_train - blup_train.mean()
delta_T_valid = blup_valid - blup_valid.mean()

def calc_C(delta_row, r_mat):
    C = 0
    n = len(delta_row)
    for i in range(n-1):
        for j in range(i+1, n):
            r_ij = r_mat.iloc[i, j]
            C += abs(r_ij) * (delta_row.iloc[i] - delta_row.iloc[j])**2
    return C

C_train = delta_T_train.apply(lambda row: calc_C(row, r_matrix), axis=1)
C_valid = delta_T_valid.apply(lambda row: calc_C(row, r_matrix), axis=1)

delta_S_train = delta_T_train.values @ omega - lambda_ * C_train.values
delta_S_valid = delta_T_valid.values @ omega - lambda_ * C_valid.values

# ============ 6. 预测精度评估 ============
rho, p_val = stats.spearmanr(delta_S_train, delta_S_valid)
rmse = np.sqrt(np.mean((delta_S_train - delta_S_valid)**2))

print("\n" + "="*50)
print(f"SLV模型预测精度（Spearman ρ）：{rho:.4f}")
print(f"P值：{p_val:.4e}")
print(f"RMSE：{rmse:.4f}")

pred_df = pd.DataFrame({'Line_ID': blup_train.index, 'delta_S': delta_S_train})
obs_df = pd.DataFrame({'Line_ID': blup_valid.index, 'delta_S_obs': delta_S_valid})
merged = pd.merge(pred_df, obs_df, on='Line_ID')

# ============ 新增：LOYO交叉验证 ============
print("\n" + "="*60)
print("LOYO 交叉验证：SLV vs Classical vs Single-trait (Isoflavone)")
print("="*60)

all_years = sorted(pheno["Year"].unique())
loyo_results = {"SLV": [], "Classical": [], "Isoflavone": []}

for valid_year in all_years:
    train_years = [y for y in all_years if y != valid_year]
    train = pheno[pheno["Year"].isin(train_years)].copy()
    valid = pheno[pheno["Year"] == valid_year].copy()

    blup_train = calc_blup(train, traits)
    blup_valid = calc_blup(valid, traits)
    common = blup_train.index.intersection(blup_valid.index)
    blup_train = blup_train.loc[common]
    blup_valid = blup_valid.loc[common]

    r_matrix = blup_train.corr()
    delta_T_train = blup_train - blup_train.mean()
    delta_T_valid = blup_valid - blup_valid.mean()

    C_train = delta_T_train.apply(lambda row: calc_C(row, r_matrix), axis=1)
    C_valid = delta_T_valid.apply(lambda row: calc_C(row, r_matrix), axis=1)

    delta_S_train = delta_T_train.values @ omega - lambda_ * C_train.values
    delta_S_valid = delta_T_valid.values @ omega - lambda_ * C_valid.values

    # Classical: λ = 0
    delta_S_train_classical = delta_T_train.values @ omega
    delta_S_valid_classical = delta_T_valid.values @ omega

    # Single-trait: Isoflavone
    iso_train = blup_train["Isoflavone"].values
    iso_valid = blup_valid["Isoflavone"].values

    rho_slv, _ = stats.spearmanr(delta_S_train, delta_S_valid)
    rho_classical, _ = stats.spearmanr(delta_S_train_classical, delta_S_valid_classical)
    rho_iso, _ = stats.spearmanr(iso_train, iso_valid)

    loyo_results["SLV"].append(rho_slv)
    loyo_results["Classical"].append(rho_classical)
    loyo_results["Isoflavone"].append(rho_iso)

    print(f"留出年份 {valid_year}: SLV ρ={rho_slv:.4f}, Classical ρ={rho_classical:.4f}, Isoflavone ρ={rho_iso:.4f}")

for key in loyo_results:
    arr = np.array(loyo_results[key])
    print(f"\n{key} LOYO mean ρ = {arr.mean():.4f} ± {arr.std(ddof=1):.4f}")

# ============ 新增：验证年 Classical 和单性状 ρ ============
print("\n" + "="*60)
print("验证年（2015）Classical 和单性状 ρ")
print("="*60)

# 使用全训练集（2006,2011,2012,2014）和验证集2015
train_full = pheno[pheno["Year"].isin([2006, 2011, 2012, 2014])].copy()
valid_full = pheno[pheno["Year"] == 2015].copy()

blup_train_full = calc_blup(train_full, traits)
blup_valid_full = calc_blup(valid_full, traits)
common_full = blup_train_full.index.intersection(blup_valid_full.index)
blup_train_full = blup_train_full.loc[common_full]
blup_valid_full = blup_valid_full.loc[common_full]

r_matrix_full = blup_train_full.corr()
delta_T_train_full = blup_train_full - blup_train_full.mean()
delta_T_valid_full = blup_valid_full - blup_valid_full.mean()

C_train_full = delta_T_train_full.apply(lambda row: calc_C(row, r_matrix_full), axis=1)
C_valid_full = delta_T_valid_full.apply(lambda row: calc_C(row, r_matrix_full), axis=1)

delta_S_train_full = delta_T_train_full.values @ omega - lambda_ * C_train_full.values
delta_S_valid_full = delta_T_valid_full.values @ omega - lambda_ * C_valid_full.values

delta_S_train_classical = delta_T_train_full.values @ omega
delta_S_valid_classical = delta_T_valid_full.values @ omega

iso_train_full = blup_train_full["Isoflavone"].values
iso_valid_full = blup_valid_full["Isoflavone"].values

rho_classical, p_classical = stats.spearmanr(delta_S_train_classical, delta_S_valid_classical)
rho_iso_valid, p_iso_valid = stats.spearmanr(iso_train_full, iso_valid_full)

print(f"Classical 验证年 ρ={rho_classical:.4f}, P={p_classical:.4e}")
print(f"单性状异黄酮 验证年 ρ={rho_iso_valid:.4f}, P={p_iso_valid:.4e}")

# ============ 7. 补画缺失的图表 ============
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 图1：SLV预测 vs 实际散点图
plt.figure(figsize=(6, 5))
plt.scatter(merged['delta_S'], merged['delta_S_obs'], alpha=0.6, color='steelblue')
plt.plot([merged['delta_S'].min(), merged['delta_S'].max()],
         [merged['delta_S'].min(), merged['delta_S'].max()],
         'r--', linewidth=2)
plt.xlabel('Predicted ΔS (Training Set)')
plt.ylabel('Observed ΔS (Validation Set)')
plt.title(f'SLV Prediction (Spearman ρ = {rho:.3f})')
plt.tight_layout()
plt.savefig("SLV_prediction_scatter.pdf")
plt.savefig("SLV_prediction_scatter.png", dpi=300)
plt.close()

# 图2：方法对比柱状图（使用实际计算值）
comparison = pd.DataFrame({
    'Method': ['SLV (ΔS)', 'Single-trait (Isoflavone)', 'Classic Breeder'],
    'Rho': [rho, rho_iso_valid, rho_classical]
})
methods = comparison['Method']
rhos = comparison['Rho']
plt.figure(figsize=(8, 5))
bars = plt.bar(methods, rhos, color=['#4C72B0', '#DD8452', '#55A868'])
plt.ylabel('Spearman ρ')
plt.title('Prediction Accuracy: SLV vs. Traditional Methods')
for bar, r in zip(bars, rhos):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
             f'{r:.4f}', ha='center', va='bottom')
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig("method_comparison_barplot.pdf")
plt.savefig("method_comparison_barplot.png", dpi=300)
plt.close()

print("图表已重新生成：SLV_prediction_scatter 和 method_comparison_barplot")
merged.to_excel("SLV_Detailed_Predictions.xlsx", index=False)
print("预测明细已保存为 SLV_Detailed_Predictions.xlsx")

# ============ 帕累托前沿分析（修正版） ============
print("\n" + "="*60)
print("帕累托前沿分析：SLV vs 单性状选择")
print("="*60)

train_pred = pd.DataFrame({
    'Line_ID': blup_train.index,
    'delta_S_train': delta_S_train
}).set_index('Line_ID')

train_iso = blup_train[['Isoflavone']].copy()
train_iso.columns = ['Isoflavone_train']

train_traits = blup_train[['Yield', 'Isoflavone']].copy()
train_traits.columns = ['Yield_train', 'Isoflavone_train']
train_selection = train_pred.join(train_traits, how='inner')

slv_top10_ids = train_selection.nlargest(10, 'delta_S_train').index.tolist()
iso_top10_ids = train_selection.nlargest(10, 'Isoflavone_train').index.tolist()

print(f"\nSLV Top10 家系：{slv_top10_ids}")
print(f"\n单性状 Top10 家系：{iso_top10_ids}")

valid_df = blup_valid.copy()
valid_df['Line_ID'] = valid_df.index

slv_valid = valid_df[valid_df['Line_ID'].isin(slv_top10_ids)].copy()
slv_valid['Group'] = 'SLV_Top10'

iso_valid = valid_df[valid_df['Line_ID'].isin(iso_top10_ids)].copy()
iso_valid['Group'] = 'SingleTrait_Top10'

comparison_df = pd.concat([slv_valid, iso_valid], ignore_index=True)

slv_yield_mean = slv_valid['Yield'].mean()
iso_yield_mean = iso_valid['Yield'].mean()
slv_iso_mean = slv_valid['Isoflavone'].mean()
iso_iso_mean = iso_valid['Isoflavone'].mean()

print(f"\n--- 验证集（2015）表现对比 ---")
print(f"SLV Top10 产量均值：{slv_yield_mean:.4f}")
print(f"单性状 Top10 产量均值：{iso_yield_mean:.4f}")
print(f"差异：{slv_yield_mean - iso_yield_mean:+.4f}")

print(f"\nSLV Top10 异黄酮均值：{slv_iso_mean:.4f}")
print(f"单性状 Top10 异黄酮均值：{iso_iso_mean:.4f}")
print(f"差异：{slv_iso_mean - iso_iso_mean:+.4f}")

from scipy import stats
t_yield, p_yield = stats.ttest_ind(slv_valid['Yield'], iso_valid['Yield'])
t_iso, p_iso = stats.ttest_ind(slv_valid['Isoflavone'], iso_valid['Isoflavone'])

print(f"\n产量差异 t检验：t={t_yield:.3f}, P={p_yield:.4f}")
print(f"异黄酮差异 t检验：t={t_iso:.3f}, P={p_iso:.4f}")

def permutation_test(group1, group2, n_perm=10000):
    observed_diff = np.mean(group1) - np.mean(group2)
    combined = np.concatenate([group1, group2])
    count = 0
    for _ in range(n_perm):
        np.random.shuffle(combined)
        perm_diff = np.mean(combined[:len(group1)]) - np.mean(combined[len(group1):])
        if abs(perm_diff) >= abs(observed_diff):
            count += 1
    return count / n_perm

p_perm_yield = permutation_test(slv_valid['Yield'].values, iso_valid['Yield'].values)
p_perm_iso = permutation_test(slv_valid['Isoflavone'].values, iso_valid['Isoflavone'].values)

print(f"\n产量差异 置换检验 P={p_perm_yield:.4f}")
print(f"异黄酮差异 置换检验 P={p_perm_iso:.4f}")

def is_pareto_efficient(costs):
    is_efficient = np.ones(costs.shape[0], dtype=bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient] < c, axis=1) | np.all(costs[is_efficient] == c, axis=1)
    return is_efficient

costs = -valid_df[['Yield', 'Isoflavone']].values
pareto_mask = is_pareto_efficient(costs)
valid_df['is_pareto'] = pareto_mask

n_pareto_total = valid_df['is_pareto'].sum()
n_slv_pareto = valid_df[valid_df['Line_ID'].isin(slv_top10_ids)]['is_pareto'].sum()
n_iso_pareto = valid_df[valid_df['Line_ID'].isin(iso_top10_ids)]['is_pareto'].sum()

print(f"\n--- 帕累托前沿统计 ---")
print(f"总家系数：{len(valid_df)}")
print(f"帕累托前沿家系数：{n_pareto_total}")
print(f"SLV Top10 中落在帕累托前沿：{n_slv_pareto}/10")
print(f"单性状 Top10 中落在帕累托前沿：{n_iso_pareto}/10")

plt.figure(figsize=(9, 7))
plt.scatter(valid_df['Yield'], valid_df['Isoflavone'],
            c='lightgray', alpha=0.5, s=40, label='All 118 lines')
pareto_points = valid_df[valid_df['is_pareto']]
plt.scatter(pareto_points['Yield'], pareto_points['Isoflavone'],
            facecolors='none', edgecolors='green', s=120, linewidth=2, label='Pareto front')
plt.scatter(slv_valid['Yield'], slv_valid['Isoflavone'],
            c='red', s=100, label='SLV Top10', zorder=5, edgecolors='black')
plt.scatter(iso_valid['Yield'], iso_valid['Isoflavone'],
            c='blue', s=100, marker='s', label='Single-trait Top10', zorder=5, edgecolors='black')

plt.xlabel('Yield BLUP (Validation Set 2015)', fontsize=12)
plt.ylabel('Isoflavone BLUP (Validation Set 2015)', fontsize=12)
plt.title('Pareto Front Analysis: SLV vs Single-trait Selection', fontsize=13)
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("pareto_front_analysis.pdf")
plt.savefig("pareto_front_analysis.png", dpi=300)
plt.close()

print("\n帕累托前沿图已保存为 pareto_front_analysis.pdf / .png")
print("="*60)

# ============ 效应量计算与敏感性分析（修改：加入P值和保存） ============
print("\n" + "="*60)
print("效应量与敏感性分析")
print("="*60)

def cohens_d(g1, g2):
    n1, n2 = len(g1), len(g2)
    var1, var2 = np.var(g1, ddof=1), np.var(g2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    return (np.mean(g1) - np.mean(g2)) / pooled_std

d_yield = cohens_d(slv_valid['Yield'].values, iso_valid['Yield'].values)
d_iso = cohens_d(slv_valid['Isoflavone'].values, iso_valid['Isoflavone'].values)
print(f"产量效应量 (Cohen's d): {d_yield:.3f}")
print(f"异黄酮效应量 (Cohen's d): {d_iso:.3f}")

# ============ 修改：Top10/20/30 敏感性分析，加入P值并保存 ============
top_results = []
for top_n in [10, 20, 30]:
    slv_top = train_selection.nlargest(top_n, 'delta_S_train').index.tolist()
    iso_top = train_selection.nlargest(top_n, 'Isoflavone_train').index.tolist()
    slv_val = valid_df[valid_df['Line_ID'].isin(slv_top)]
    iso_val = valid_df[valid_df['Line_ID'].isin(iso_top)]

    # 置换检验
    p_perm_yield = permutation_test(slv_val['Yield'].values, iso_val['Yield'].values)
    p_perm_iso = permutation_test(slv_val['Isoflavone'].values, iso_val['Isoflavone'].values)

    top_results.append({
        'Top_n': top_n,
        'SLV_Yield_mean': slv_val['Yield'].mean(),
        'SLV_Iso_mean': slv_val['Isoflavone'].mean(),
        'Iso_Yield_mean': iso_val['Yield'].mean(),
        'Iso_Iso_mean': iso_val['Isoflavone'].mean(),
        'Yield_diff': slv_val['Yield'].mean() - iso_val['Yield'].mean(),
        'Iso_diff': slv_val['Isoflavone'].mean() - iso_val['Isoflavone'].mean(),
        'P_yield': p_perm_yield,
        'P_iso': p_perm_iso
    })

    print(f"\nTop{top_n}：")
    print(f"  SLV 产量均值 = {slv_val['Yield'].mean():.4f}，异黄酮均值 = {slv_val['Isoflavone'].mean():.4f}")
    print(f"  单性状 产量均值 = {iso_val['Yield'].mean():.4f}，异黄酮均值 = {iso_val['Isoflavone'].mean():.4f}")
    print(f"  产量差异 = {slv_val['Yield'].mean() - iso_val['Yield'].mean():+.4f}，P = {p_perm_yield:.4f}")
    print(f"  异黄酮差异 = {slv_val['Isoflavone'].mean() - iso_val['Isoflavone'].mean():+.4f}，P = {p_perm_iso:.4f}")

top_df = pd.DataFrame(top_results)
top_df.to_excel("Top_selection_results.xlsx", index=False)
print("\nTop10/20/30 结果已保存到 Top_selection_results.xlsx")
print(top_df)

# ============ 补充分析1：系统整合度 S ============
print("\n" + "="*60)
print("系统整合度 S 计算")
print("="*60)

def calc_entropy(series, bins=10):
    hist, _ = np.histogram(series, bins=bins, density=True)
    hist = hist[hist > 0] * np.diff(_)[0]
    hist = hist / hist.sum()
    return -np.sum(hist * np.log2(hist))

def calc_joint_entropy(df, cols, bins=10):
    hist, edges = np.histogramdd(df[cols].values, bins=bins, density=True)
    volumes = np.prod([np.diff(e)[0] for e in edges])
    probs = hist.flatten() * volumes
    probs = probs[probs > 0]
    probs = probs / probs.sum()
    return -np.sum(probs * np.log2(probs))

trait_cols = ['Yield', 'Isoflavone', 'Seed_Weight']
H_joint = calc_joint_entropy(blup_valid, trait_cols)
H_marginal = sum(calc_entropy(blup_valid[col]) for col in trait_cols)
S = 1 - H_joint / H_marginal

print(f"联合熵 H_joint = {H_joint:.4f}")
print(f"边际熵之和 H_max = {H_marginal:.4f}")
print(f"系统整合度 S = {S:.4f}")

blup_valid_std = blup_valid[trait_cols].std(axis=1)
blup_valid_mean = blup_valid[trait_cols].mean(axis=1)
cv = blup_valid_std / blup_valid_mean.abs()
blup_valid['S_value'] = S
blup_valid['CV'] = cv

# ============ 补充分析2：道地适应性 A(g,e) ============
print("\n" + "="*60)
print("道地适应性 A(g,e) 计算")
print("="*60)

train['composite'] = train['Yield'] * 0.35 + train['Isoflavone'] * 0.45 + train['Seed_Weight'] * 0.20

f_g = train.groupby('Line_ID')['composite'].mean()
f_e = train.groupby('Location')['composite'].mean()
grand_mean = train['composite'].mean()

results = []
for line in f_g.index:
    for loc in f_e.index:
        sub = train[(train['Line_ID'] == line) & (train['Location'] == loc)]
        if len(sub) == 0:
            continue
        obs_mean = sub['composite'].mean()
        I_ge = obs_mean - f_g[line] - f_e[loc] + grand_mean
        A = (f_g[line] * f_e[loc] * I_ge) / (f_g[line] + f_e[loc]) if (f_g[line] + f_e[loc]) != 0 else 0
        results.append({'Line_ID': line, 'Location': loc, 'A_ge': A})

A_df = pd.DataFrame(results)
A_pivot = A_df.pivot(index='Line_ID', columns='Location', values='A_ge')
print("A(g,e) 计算完成，各地点适应性均值：")
print(A_pivot.mean())

A_pivot.to_excel("A_ge_results.xlsx")

# ============ 补充分析3：产量导向Top10对照组 ============
print("\n" + "="*60)
print("产量导向Top10对照组分析")
print("="*60)

yield_top10_ids = train_selection.nlargest(10, 'Yield_train').index.tolist()
yield_valid = valid_df[valid_df['Line_ID'].isin(yield_top10_ids)]

print(f"产量导向Top10 家系：{yield_top10_ids}")
print(f"\n产量导向Top10 在验证集的产量均值：{yield_valid['Yield'].mean():.4f}")
print(f"产量导向Top10 在验证集的异黄酮均值：{yield_valid['Isoflavone'].mean():.4f}")

print("\n" + "="*60)
print("三组对比汇总（验证集2015）")
print("="*60)
print(f"{'组别':<20} {'产量均值':<12} {'异黄酮均值':<12}")
print(f"{'SLV Top10':<20} {slv_valid['Yield'].mean():<12.4f} {slv_valid['Isoflavone'].mean():<12.4f}")
print(f"{'异黄酮导向Top10':<20} {iso_valid['Yield'].mean():<12.4f} {iso_valid['Isoflavone'].mean():<12.4f}")
print(f"{'产量导向Top10':<20} {yield_valid['Yield'].mean():<12.4f} {yield_valid['Isoflavone'].mean():<12.4f}")

plt.figure(figsize=(8, 6))
plt.scatter(slv_valid['Yield'], slv_valid['Isoflavone'], c='red', s=100, label='SLV Top10', zorder=5)
plt.scatter(iso_valid['Yield'], iso_valid['Isoflavone'], c='blue', s=100, marker='s', label='Isoflavone Top10', zorder=5)
plt.scatter(yield_valid['Yield'], yield_valid['Isoflavone'], c='green', s=100, marker='^', label='Yield Top10', zorder=5)
plt.xlabel('Yield BLUP (Validation Set 2015)')
plt.ylabel('Isoflavone BLUP (Validation Set 2015)')
plt.title('Three-way Comparison: SLV vs Single-trait Selection')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("three_way_comparison.pdf")
plt.savefig("three_way_comparison.png", dpi=300)
plt.close()

print("\n三组对比图已保存为 three_way_comparison.pdf / .png")
print("="*60)
print("\n分析完成！")