# ============ 补充诊断：单性状预测精度 ============
print("\n" + "="*50)
print("单性状预测精度诊断：")
print("="*50)

for i, tr in enumerate(traits):
    rho_single, _ = stats.spearmanr(blup_train[tr], blup_valid[tr])
    print(f"  {tr} 单性状预测 ρ = {rho_single:.4f}")

# 查看异黄酮与产量在验证集的联合分布
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for i, tr in enumerate(traits):
    axes[i].scatter(blup_train[tr], blup_valid[tr], alpha=0.6)
    axes[i].set_xlabel(f'{tr} (训练集)')
    axes[i].set_ylabel(f'{tr} (验证集)')
    axes[i].set_title(f'{tr} 预测 vs 实际')
plt.tight_layout()
plt.savefig("single_trait_scatter.png")
print("散点图已保存为 single_trait_scatter.png")