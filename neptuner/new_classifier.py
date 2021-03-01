import sys
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn import metrics
from sklearn.model_selection import train_test_split
import pickle

###-MAIN FUNCTION-###
## USAGE: python3 new_classifier.py ndrname exprname


def new_classifier(ndrname, exprname):
    # loading tables
    ndrt = pd.read_csv(ndrname)
    labels = pd.read_csv(exprname)

    # preparation in the form features|label
    wtexpr = np.array([])
    for i in np.arange(ndrt.shape[0]):
        row = int(np.where(labels["Names"] == ndrt.loc[i, "Name"])[0])
        wtexpr = np.append(wtexpr, labels.loc[row, "WTexp"])

    ndrt["wtexpr"] = wtexpr

    # setup reference table
    ndrt = ndrt[
        ndrt[["plus1length", "minus1length", "plus1height", "minus1height", "ndr"]].notnull().all(1)]
    ndrt = ndrt[ndrt.notnull().all(1)]

    # svm parameters
    X = ndrt[["plus1length", "minus1length", "plus1height", "minus1height", "ndr"]]
    y = ndrt["wtexpr"]

    rs = 123
    weights = {0: (len(y) - sum(y)), 1: sum(y)}

    # training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        stratify=y,
                                                        test_size=0.1,
                                                        random_state=123)

    # GRID SEARCH
    from sklearn.model_selection import GridSearchCV
    param_grid = {'C': [0.1, 1, 10, 100], 'gamma': [1, 0.1, 0.01, 0.001], 'kernel': ['rbf']}
    grid = GridSearchCV(svm.SVC(), param_grid, refit=True, verbose=2)
    grid.fit(X_train, y_train)
    print(grid.best_estimator_)
    best_y_pred = grid.predict(X_test)
    print(metrics.precision_score(y_test, best_y_pred))  # Output

    filename = 'finalized_model.sav'
    pickle.dump(grid, open(filename, 'wb'))

    return


###-END FUNCTIONS-###


ndrtable = sys.argv[1]
expr_labels = sys.argv[2]

new_classifier(ndrtable, expr_labels)
