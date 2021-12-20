import numpy as np

def find_distinct_values(df,attributes_lst):
    distinct_values = {}
    for attribute in attributes_lst:
        distinct_values[attribute] = set(df[attribute])

    return distinct_values
	
def calculate_Laplace_smoothed_marginals(df,distinct_values):
    #initialize the marginals
    marginals = {}
    for attribute_type in distinct_values:
        for attribute in distinct_values[attribute_type]:
            marginals[attribute] = np.array([1,1])  #initializing to [1,1] implements Laplace smoothing
            
    for attribute_type in distinct_values:
        for attribute, authenticity in zip(df[attribute_type],df['authenticity']):
            marginals[attribute][authenticity] += 1
            
    return marginals    
	
def characterize_data_set(df):
    fake_label           = 0
    true_label           = 1
    summary              = {}
    authenticity_data    = df['authenticity']
    fake_counts          = len(np.where(authenticity_data==fake_label)[0])
    true_counts          = len(np.where(authenticity_data==true_label)[0])
    #counts               = authenticity_data.value_counts()
	#if len(counts) < 2:
	#    print(counts)
	#    #data of all one type - let's find out which one
	#	authenticity_val = df['authenticity'][0]
	#	if authenticity_val == 0:
	#	    counts = [counts,0]
	#	else:
	#	    counts = [0,counts]
		
    
    summary['num samples'] = fake_counts + true_counts
    summary['num fake']    = fake_counts
    summary['num true']    = true_counts
    
    return summary
	
def calculate_evidence(distinct_values,smoothed_marginals,summary,sample_values):
    fake_label           = 0
    true_label           = 1
    num_attributes       = len(distinct_values)
    
    prob_fake            = summary['num fake']/summary['num samples']
    prob_true            = summary['num true']/summary['num samples']
    smoothed_num_fake    = summary['num fake'] + num_attributes
    smoothed_num_true    = summary['num true'] + num_attributes
    
    sample_evidence_fake = 1
    for attribute in sample_values:
        sample_evidence_fake *= smoothed_marginals[attribute][fake_label]/smoothed_num_fake
    sample_evidence_fake *= prob_fake
    
    sample_evidence_true = 1
    for attribute in sample_values:
        sample_evidence_true *= smoothed_marginals[attribute][true_label]/smoothed_num_true
    sample_evidence_true *= prob_true
    
    normalization = sample_evidence_fake + sample_evidence_true
    
    return sample_evidence_fake/normalization, sample_evidence_true/normalization
	
