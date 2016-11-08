import matplotlib.pyplot as plt
from sklearn import svm, metrics

def roc_curve(actual_values, predicted_values):
	fpr, tpr, thresholds = metrics.roc_curve(actual_values, predicted_values, pos_label=1)
	roc_auc = metrics.roc_auc_score(actual_values, predicted_values, average='macro', sample_weight=None)
	plt.title('Receiver Operating Characteristic')
	plt.plot(fpr, tpr, 'b', label='AUC = %0.2f'% (roc_auc))
	plt.legend(loc='lower right')
	plt.plot([0,1],[0,1],'r--')
	plt.xlim([-0.1,1.2])
	plt.ylim([-0.1,1.2])
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	return plt
