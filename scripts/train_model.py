import pandas as pd
import joblib
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

# Load merged data
df = pd.read_csv(snakemake.input[0])

# Split X and Y
X = df.drop(['molecule_chembl_id', 'pIC50'], axis=1)
Y = df['pIC50']

# Drop NaNs
data = pd.concat([X, Y], axis=1).dropna()
X = data.drop('pIC50', axis=1)
Y = data['pIC50']

# Train-test split
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# Train model
model = RandomForestRegressor(n_estimators=100, random_state=100)
model.fit(X_train, Y_train)

# Save model
joblib.dump(model, snakemake.output[0])