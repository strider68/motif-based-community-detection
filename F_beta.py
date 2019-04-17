#Compute F-score
def precision(y_true,y_pred):
    tp=0.0
    fp=0.0
    n1=len(y_true)
    for i in range(n1):
        for j in range(i+1,n1):
            if y_true[i]==y_true[j] and y_pred[i]==y_pred[j]:
                tp+=1
            if y_true[i]!=y_true[j] and y_pred[i]==y_pred[j]:
                fp+=1
    if tp==0 and fp==0:
        return 0
    else:
        return tp/(tp+fp)

def recall(y_true,y_pred):
    tp=0.0
    fn=0.0
    n1=len(y_true)
    for i in range(n1):
        for j in range(i+1,n1):
            if y_true[i]==y_true[j] and y_pred[i]==y_pred[j]:
                tp+=1
            if y_true[i]==y_true[j] and y_pred[i]!=y_pred[j]:
                fn+=1
    if tp==0 and fn==0:
        return 0
    else:
        return tp/(tp+fn)


##--computation of F_beta,
def F_score(y_true,y_pred,beta):
    P=precision(y_true,y_pred)
    R=recall(y_true,y_pred)
    if P==0 and R==0:
        return 0
    else:
        return ((beta*beta+1)*P*R)/(beta*beta*P+R)
