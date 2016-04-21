#################################################################################
# Copyright [2016] Novartis Institutes for BioMedical Research Inc.             #
# Licensed under the Apache License, Version 2.0 (the "License");               #
# you may not use this file except in compliance with the License.              #
# You may obtain a copy of the License at                                       #
# http://www.apache.org/licenses/LICENSE-2.0                                    #
# Unless required by applicable law or agreed to in writing, software           #
# distributed under the License is distributed on an "AS IS" BASIS,             #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.      #
# See the License for the specific language governing permissions and           #
# limitations under the License.                                                #
#################################################################################

inputdata<-read.csv("example_inputdata_toolscore_calculation.csv", header=TRUE)
#drugbank_names: compound identifier http://www.drugbank.ca/
#inchi_key: compound identifier
#gene_id/symbol: target identifier http://www.ncbi.nlm.nih.gov/gene
#homologene_group_id: target identifier across species http://www.ncbi.nlm.nih.gov/homologene
#micromolar_activity: dose-dependent measurements such as IC50/EC50 between compound and target
#origin: source where data is integrated

example_query_compound<-"Nilotinib"
example_query_target<-"ABL1"
#use Nilotinib and ABL1 as tool_score calculation example

calc_tool_score<-function(inputdata, example_query_compound, example_query_target) 
{
  example_query_homologene_group_id<-inputdata[inputdata$symbol==example_query_target,]$homologene_group_id[1]
  example_subset<-inputdata[inputdata$drugbank_names==example_query_compound,]
  
  ontarget<-example_subset[example_subset$homologene_group_id == example_query_homologene_group_id,]
  offtarget<-example_subset[example_subset$homologene_group_id != example_query_homologene_group_id,]
  
  drugbank_active<-ontarget[ontarget$origin=='drugbank',]
  drugbank_strength<-7 #weight 7 for approved drugs, 2 for approved nutraceuticals, 1 else
  
  chembl_active<-ontarget[ontarget$origin=='chembl',]
  chembl_active_data<-as.numeric(as.character(chembl_active$micromolar_activity))
  chembl_active_low_nM<-chembl_active_data[chembl_active_data<=0.1]
  if (length(chembl_active_low_nM) > 1) {
    chembl_strength<-3
  }else if (length(chembl_active_data[chembl_active_data<=1]) > 4) {
    chembl_strength<-2
  }else {
    chembl_strength<-1
  }
  
  strength<-max(chembl_strength, drugbank_strength) + 1 # bonus for multiple sources
  
  ontarget_IC50<-as.numeric(as.character(ontarget[!is.na(ontarget$micromolar_activity),]$micromolar_activity))
  ontarget_IC50_N<-length(ontarget_IC50)
  ontarget_IC50_Q1<-quantile(ontarget_IC50, probs=0.25, na.rm=TRUE)
  offtarget_IC50<-as.numeric(as.character(offtarget[!is.na(offtarget$micromolar_activity),]$micromolar_activity))
  offtarget_IC50_N<-length(offtarget_IC50)
  offtarget_IC50_Q1<-quantile(offtarget_IC50, probs=0.25, na.rm=TRUE)
  Q1_IC50_diff<-log10(offtarget_IC50_Q1)-log10(ontarget_IC50_Q1);
  if (length(ontarget_IC50)==0 || length(offtarget_IC50)==0) {
    wilcox_pval<-10;
  } else {
    w<-wilcox.test(ontarget_IC50, offtarget_IC50, alternative='less');
    if (w$p.value>=1E-32){
      wilcox_pval<-w$p.value;
    } else {
      wilcox_pval<-1E-32;
    }
  }
  
  IC50<-c(ontarget_IC50, offtarget_IC50)
  IC50_Q1<-quantile(IC50, probs=0.25, na.rm=TRUE)
  investigation_bias_Q1<-length(ontarget_IC50[ontarget_IC50<IC50_Q1])/length(IC50[IC50<IC50_Q1])
  
  selectivity<-(Q1_IC50_diff/10 + investigation_bias_Q1 - log10(wilcox_pval)/32)/3
  tool_score<- strength * selectivity
  names(tool_score)<-'tool score'
  return (tool_score)
}

calc_tool_score(inputdata, example_query_compound, example_query_target)
