import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('q3_data.csv')

# Compute and print y-axis standard deviation
y_std = df['y'].std()
y_mean = df['y'].mean()
print(f"Y-axis mean : {y_mean:.6f} & standard deviation: {y_std:.6f}")

# Optional: preview the data
print(df.head())

# Use seaborn for plotting
sns.set(style="whitegrid")
ax = sns.scatterplot(data=df, x='x', y='y')
ax.set_title('Scatterplot of x vs y')
ax.set_xlabel('x')
ax.set_ylabel('y')

plt.tight_layout()
plt.show()