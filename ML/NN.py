import numpy as np
import pandas as pd
from sklearn import preprocessing, model_selection
import torch
import torch.nn.functional as F
import numpy as np
from torch.utils import data
from torch import nn
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelBinarizer
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.model_selection import StratifiedKFold

fig, c_ax = plt.subplots(1,1, figsize = (12, 8))
#target = ['CL', 'ME', 'PN']
target = ['IDHmut-codel', 'IDHmut-non-codel', 'IDHwt']
#target = ['less than 1', '1-3', 'more than 3']
EPOCH = 30
batch_size = 30
input_size = 957
#1507
#646
#957
criterion = nn.CrossEntropyLoss()


class MyDataSet(data.Dataset):
    def __init__(self, training, labels):
        super(MyDataSet, self).__init__()
        self.train = training
        self.labels = labels

    def __len__(self):
        # get number of instances and labels
        return len(self.train)

    def __getitem__(self, idx):
        # return instance and it's label
        return self.train[idx], self.labels[idx]


class Model(torch.nn.Module):
    def __init__(self, size):
        super(Model, self).__init__()
        self.batch1 = torch.nn.BatchNorm1d(100)
        self.fc0 = torch.nn.Linear(input_size, 100)
        self.batch2 = torch.nn.BatchNorm1d(50)
        self.fc1 = torch.nn.Linear(100, 50)
        self.fc2 = torch.nn.Linear(50, 3)

    def forward(self, x):
        x = x.view(-1, input_size)
        x = F.relu(self.batch1(self.fc0(x)))
        x = F.relu(self.batch2(self.fc1(x)))
        x = self.fc2(x)
        return x


def train(model, dataloader):
    model.train()
    model.double()
    train_loader = torch.utils.data.DataLoader(dataloader, batch_size=17, shuffle=True)
    for i in range(EPOCH):
        tot_loss= 0.
        for batch_idx, (data_x, labels) in enumerate(train_loader):
            optimizer.zero_grad()
            data_x = data_x.double()
            output = model(data_x)
            # labels = torch.unsqueeze(labels, dim=1)
            loss = criterion (output, labels.long())
            tot_loss +=loss
            #print(criterion(output, labels.long()).item())
            loss.backward()
            optimizer.step()
        print(tot_loss)

def test(model, test_loader):
    model.eval()
    model.double()
    test_loss = 0.
    correct = 0
    preds=[]
    all_labels=[]
    test_loader = torch.utils.data.DataLoader(test_loader, batch_size=1, shuffle=True)
    for batch_idx, (data_x, labels) in enumerate(test_loader):
        data_x = data_x.reshape(-1)
        output = model(data_x.double())
        test_loss += criterion(output, labels.long()).item()
        pred = output.max(1,keepdim=True)[1]
        correct += pred.eq(labels.view_as(pred)).sum().item()
        preds.append(output[0].tolist())
        all_labels.append(labels[0].tolist())
    test_loss /= len(test_loader)
    correct = 100. * correct / len(test_loader)
    return correct,all_labels,preds

def read_data():
    df = pd.read_csv('TCGA_DATA- IDH astro.csv')
    df = df.rename(columns={'Unnamed: 0': 'gene'})
    df = df.iloc[:, 1:]
    df.columns.values[0] = 'label'
    #df = df[df['label'] != 'NE']
    le = preprocessing.LabelEncoder()
    df['label'] = le.fit_transform(df['label'])
    df = df.to_numpy()
    np.random.shuffle(df)
    df[:, 2:] =preprocessing.minmax_scale(df[:, 2:], feature_range=(0, 1), axis=0, copy=True)
    return df



df = read_data()
X = df[:,1:]
y = df[:,0]
k = 4
tot_correct=0.
kf = KFold(n_splits=k, random_state=None)
for train_index , test_index in kf.split(X):
    X_train , X_test = X[train_index,:],X[test_index,:]
    y_train , y_test = y[train_index] , y[test_index]
    model = Model(input_size)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
    train_data = MyDataSet(X_train, y_train)
    train(model, train_data)
    test_data = MyDataSet(X_test, y_test)
    correct,labels,preds= test(model, test_data)
    tot_correct += correct
print(tot_correct/k)
labels = np.array(labels)
labels = labels.astype(int)
labels=np.eye(3)[labels]
#labels= labels.tolist()
preds = np.array(preds)

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(3):
    fpr[i], tpr[i], _ = roc_curve(labels[:, i], preds[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])


for (idx, c_label) in enumerate(target):
    fpr, tpr, thresholds = roc_curve(labels[:,idx].astype(int), preds[:,idx])
    plt.plot(fpr, tpr, label = '%s (AUC:%0.2f)'  % (c_label, auc(fpr, tpr)))
plt.plot(fpr, fpr, 'b-', label = 'Random Guessing')


plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
#plt.title('Some extension of Receiver operating characteristic to multi-class')
plt.legend(loc="lower right")
plt.show()

