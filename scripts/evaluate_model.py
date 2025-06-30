import pandas as pd
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split

# Load data and model
df = pd.read_csv(snakemake.input[0])
model = joblib.load(snakemake.input[1])

# Prepare features and target
X = df.drop(['molecule_chembl_id', 'pIC50'], axis=1)
Y = df['pIC50']

# Remove missing values
data = pd.concat([X, Y], axis=1).dropna()
X = data.drop('pIC50', axis=1)
Y = data['pIC50']

# Split the data
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# Predict
Y_pred = model.predict(X_test)

# Plot
sns.set(style='white')
plt.figure(figsize=(5,5))
ax = sns.regplot(x=Y_test, y=Y_pred, scatter_kws={'alpha':0.4})
ax.set_xlabel('Experimental pIC50')
ax.set_ylabel('Predicted pIC50')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
plt.tight_layout()
plt.savefig(snakemake.output[0])
