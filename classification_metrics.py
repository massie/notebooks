from sklearn.metrics import confusion_matrix, accuracy_score, classification_report

def classification_metrics(actual, predicted, target_names):
	print "Accuracy was %.2f%%\n" % (100 * accuracy_score(actual, predicted))
	print classification_report(actual, predicted, target_names=target_names)
	cm = confusion_matrix(actual, predicted, labels=range(len(target_names)))
	print "Confusion Matrix: cols = predictions, rows = actual\n"
	row_format ="{:>15}" * (len(target_names)+1)
	print row_format.format("", *target_names)
	for target_name, row in zip(target_names, cm):
		print row_format.format(target_name, *row)
