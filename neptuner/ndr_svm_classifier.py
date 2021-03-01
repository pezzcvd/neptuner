import sys
import pandas as pd
import numpy as np
import pickle


###-MAIN FUNCTION-###
## USAGE: python3 ndr_svm_classifier.py ndrname exprname modelname


def use_svm(ndrname, exprname, modelname):
    outname = ndrname.split("/")[-1]
    outname = outname.split(".")[0]

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

    loaded_model = pickle.load(open(modelname, 'rb'))

    y_calc = loaded_model.predict(X)

    result = loaded_model.score(X, y)
    print(result)
    ndrt = ndrt.drop("wtexpr", 1)
    ndrt["NDRpattern"] = y_calc
    ndrt.to_csv(outname + "_svm.csv", index=False)

    return


###-END FUNCTIONS-###

ndrtable = sys.argv[1]
expr_labels = sys.argv[2]
svm_model = sys.argv[3]

use_svm(ndrtable, expr_labels, svm_model)
