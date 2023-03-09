from sklearn.ensemble import RandomForestClassifier
import numpy  as np
from matplotlib import pyplot as plt
from sklearn.metrics import accuracy_score, roc_curve, auc, roc_auc_score
from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn import preprocessing
from sklearn.preprocessing import LabelBinarizer, label_binarize
from sklearn.metrics import accuracy_score

fig, c_ax = plt.subplots(1,1, figsize = (12, 8))
#target = ['CL', 'ME', 'PN']
#arget = ['IDHmut-codel', 'IDHmut-non-codel', 'IDHwt']
#target = ['less than 1', '1-3', 'more than 3']
#target = ['less than 1', 'more than 3']


def read_data():
    df = pd.read_csv('TCGA_DATA- LGG old survival Heiland.csv')
    df = df.rename(columns={'Unnamed: 0': 'gene'})
    df = df.iloc[:, 1:]
    df.columns.values[0] = 'label'
    df = df[df['label'] != '1-3']
    le = preprocessing.LabelEncoder()
    df['label'] = le.fit_transform(df['label'])
    target = le.classes_
    df = df.to_numpy()
    np.random.shuffle(df)
    df[:, 2:] =preprocessing.minmax_scale(df[:, 2:], feature_range=(0, 1), axis=0, copy=True)
    return df, target

def random_forest_classifier(x_train, y_train):
    clf = RandomForestClassifier()
    clf.fit(x_train, y_train)

    return clf



df, target = read_data()
train, test = train_test_split(df, test_size=0.2, shuffle=True)
x_train = train[:,1:]
y_train = train[:,0]
x_test = test[:,1:]
y_test = test[:,0]

# Create random forest classifier instance
trained_model = random_forest_classifier(x_train, y_train)
predictions = trained_model.predict(x_test)
y_score = trained_model.predict_proba(x_test)
accuracy = accuracy_score(y_test, predictions)
print(accuracy)
#multiclass_roc_auc_score(y_test, y_score)

#Binarize the output
y_test_bin = label_binarize(y_test, classes=[0,1,2])
n_classes = y_test_bin.shape[1]

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(len(target)):
  fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_score[:, i])
  plt.plot(fpr[i], tpr[i], label='%s (AUC:%0.2f)' % (target[i], auc(fpr[i], tpr[i])))

plt.legend()
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
#plt.title('Receiver Operating Characteristic Curves')
plt.show()