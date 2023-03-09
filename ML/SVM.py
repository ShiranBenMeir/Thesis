import random
import numpy  as np
from matplotlib import pyplot as plt
from sklearn.metrics import accuracy_score, roc_curve, auc, roc_auc_score
from sklearn.model_selection import train_test_split
from PIL import Image
import pandas as pd
from sklearn import preprocessing
from sklearn.preprocessing import LabelBinarizer
from sklearn.svm import SVC
from sklearn import svm, datasets
from sklearn.preprocessing import LabelBinarizer, label_binarize
import sklearn.model_selection as model_selection
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score

#target= ['ME', 'PN', 'CL']
#target = ['less than 1', '1-3', 'more than 3']
#target = ['IDHmut-codel', 'IDHmut-non-codel', 'IDHwt']


#fig, c_ax = plt.subplots(1,1, figsize = (12  , 8))


def read_data():
    df = pd.read_csv('TCGA_DATA- LGG old survival Heiland.csv')
    df = df.rename(columns={'Unnamed: 0': 'gene'})
    df = df.iloc[:, 1:]
    df.columns.values[0] = 'label'
    df = df[df['label'] != '1-3']
    le = preprocessing.LabelEncoder()
    df['label'] = le.fit_transform(df['label'])
    target=le.classes_
    df = df.to_numpy()
    np.random.shuffle(df)
    df[:, 2:] =preprocessing.minmax_scale(df[:, 2:], feature_range=(0, 1), axis=0, copy=True)
    return df, target

def ROC_generator(y_test, y_proba, target):
    y_test_bin = label_binarize(y_test, classes=[0, 1, 2])
    n_classes = y_test_bin.shape[1]

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(2):
        fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_proba[:, i])
        plt.plot(fpr[i], tpr[i], label='%s (AUC:%0.2f)' % (target[i], auc(fpr[i], tpr[i])))

    plt.legend()
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # plt.title('Receiver Operating Characteristic Curves')
    plt.show()


def accuracy_calculation(y_test,poly_pred, rbf_pred, sigmoid_pred, linear_pred):
    # calculate accuracy
    poly_accuracy = accuracy_score(y_test, poly_pred)
    rbf_accuracy = accuracy_score(y_test, rbf_pred)
    sigmoid_accuracy = accuracy_score(y_test, sigmoid_pred)
    linear_accuracy = accuracy_score(y_test, linear_pred)
    print('polynomial accuracy : ', "%.2f" % (poly_accuracy * 100))
    print('radial accuracy : ', "%.2f" % (rbf_accuracy * 100))
    print('sigmoid accuracy : ', "%.2f" % (sigmoid_accuracy * 100))
    print('linear accuracy : ', "%.2f" % (linear_accuracy * 100))

    return poly_accuracy, rbf_accuracy, sigmoid_accuracy, linear_accuracy

def SVM_model(x_train, y_train, x_test):
    rbf_model = svm.SVC(kernel='rbf',probability=True).fit(x_train, y_train)
    poly_model = svm.SVC(kernel='poly', probability=True).fit(x_train, y_train)
    sigmoid_model = svm.SVC(kernel='sigmoid',probability=True).fit(x_train, y_train)
    linear_model = svm.SVC(kernel='linear',probability=True).fit(x_train, y_train)
    #use the model for prediction
    poly_pred = poly_model.predict(x_test)
    rbf_pred = rbf_model.predict(x_test)
    sigmoid_pred = sigmoid_model.predict(x_test)
    linear_pred = linear_model.predict(x_test)

    #y_praba = poly_model.predict_proba(x_test)
    #y_praba = rbf.predict_proba(x_test)
    #y_praba = sigmoid.predict_proba(x_test)
    #y_praba = linear.predict_proba(x_test)

    return poly_pred, rbf_pred, sigmoid_pred, linear_pred, rbf_model, poly_model, sigmoid_model, linear_model


df, target = read_data()
train, test = train_test_split(df, test_size=0.2, shuffle=True)
x_train = train[:,1:]
y_train = train[:,0]
x_test = test[:,1:]
y_test = test[:,0]
poly_pred, rbf_pred, sigmoid_pred, linear_pred, rbf_model, poly_model, sigmoid_model, linear_model =\
    SVM_model(x_train, y_train, x_test)

poly_accuracy, rbf_accuracy, sigmoid_accuracy, linear_accuracy = \
    accuracy_calculation(y_test,poly_pred, rbf_pred, sigmoid_pred, linear_pred)
accuracy_lst=[poly_accuracy, rbf_accuracy, sigmoid_accuracy, linear_accuracy]
max_accuracy_idx = accuracy_lst.index(max(accuracy_lst))

if max_accuracy_idx == 0:
    best_model= poly_model
if max_accuracy_idx == 1:
    best_model= rbf_model
if max_accuracy_idx == 2:
    best_model= sigmoid_model
if max_accuracy_idx == 3:
    best_model= linear_model

y_proba = best_model.predict_proba(x_test)
ROC_generator(y_test, y_proba, target)








