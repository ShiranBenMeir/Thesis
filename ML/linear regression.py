# Multiple Linear Regression

# Importing the libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pandas as pd
import numpy as np

# Importing the dataset
df = pd.read_csv('TCGA_DATA- old normalized numeric survival Heiland.csv')
df = df.rename(columns={'Unnamed: 0': 'gene'})
df = df.iloc[:, 1:]
df.columns.values[0] = 'label'

# Min-Max Normalization
df2 = df.drop('label', axis=1)
df_norm = (df2-df2.min())/(df2.max()-df2.min())
df_norm = pd.concat((df.label, df_norm), 1)
df = df_norm

# shuffle all rows in df
df = df.sample(frac=1).reset_index()
df= df.drop(df.columns[np.isnan(df).any()], axis=1)
X = df.iloc[:, 1:]
y = df.iloc[:, 0]


# Splitting the dataset into the Training set and Test set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)

# Fitting Multiple Linear Regression to the Training set
from sklearn.linear_model import LinearRegression
regressor = LinearRegression()
regressor.fit(X_train, y_train)

# Predicting the Test set results
y_pred = regressor.predict(X_test)

from sklearn.metrics import r2_score
score=r2_score(y_test,y_pred)
print(score)

from sklearn.metrics import mean_squared_error

mse = mean_squared_error(y_pred, y_test)
print("Mean Squarred Error is :" ,mse*100)

# plot
originalStability = 1
plt.figure()
#x = x.view(-1, input_size)
plt.scatter(y_test, y_pred)
plt.xlabel("real survival")
plt.ylabel("predictional survival")
plt.title("linear regression")
z = np.polyfit(y_test, y_pred, 1)
p = np.poly1d(z)
plt.plot(y_test, p(y_test), "r--")
# R = r2_score(preds, p(stabs))
# plt.text(4.5, 5.5, 'R^2=%0.3f' % R, fontdict={'fontsize': 17})
plt.savefig("linear regression")
plt.show()