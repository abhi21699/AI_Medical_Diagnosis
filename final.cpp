#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <bits/stdc++.h>
#include <random>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;
using namespace chrono;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<double> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name,int n,vector<string> vals)
	{
		Node_Name=name;
		nvalues=n;
		values=vals;
	}
	string get_name()
	{
		return Node_Name;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<string> get_Parents()
	{
		return Parents;
	}
	vector<double> get_CPT()
	{
		return CPT;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
	void set_CPT(vector<double> new_CPT)
	{
		CPT.clear();
		CPT=new_CPT;
	}
    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents=Parent_Nodes;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index )
    {
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }

    // A function to print the current node of the bayes net
    void printNode() {
        cout << "Variable Name: " << Node_Name << endl;
        cout << "Number of Categories: " << nvalues << endl;
        cout << "Categories: ";
        for (const string& value : values) {
            cout << value << " ";
        }
        cout << endl;

        cout << "Parents: ";
        for (const string& parent: Parents) {
            cout << parent << " ";
        }
        cout << endl;

        cout << "Children: ";
        for (const int child: Children) {
            cout << child << " ";
        }
        cout << endl;

        cout << "Conditional Probability Table (CPT): ";
        for (const double prob: CPT) {
            cout << prob << " ";
        }
        cout << endl;

        cout << "Size of CPT: " << CPT.size() << endl;
    }

    // A function to print the cpts for all the variables
    void printCPT() {
        cout << Node_Name << " ";
        for (const double prob : CPT) {
            cout << prob << " ";
        }
        cout << endl;
    }

    // The following function will return the index of the value of a node.
    int valueIndex(const string& value) {
        int index = 0;
        for (const string& curr_value : values) {
            if (value == curr_value) {
                break;
            }
            index++;
        }
        return index;
    }

};

 // The whole network represted as a list of nodes
class network{

	list <Graph_Node> Pres_Graph;
    // The following vector will store the keys to find the probabilities in the CPT
    vector<vector<int>> cpt_indexes;
    vector<vector<int>> parent_num_values;

public:
	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}
    
    
	int netSize()
	{
		return Pres_Graph.size();
	}
    // get the index of node with a given name
    int get_index(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return count;
            count++;
        }
        return -1;
    }
// get the node at nth index
    list<Graph_Node>::iterator get_nth_node(int n)
    {
       list<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(count==n)
                return listIt;
            count++;
        }
        return listIt; 
    }
    //get the iterator of a node with a given name
    list<Graph_Node>::iterator search_node(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return listIt;
        }
    
            cout<<"node not found\n";
        return listIt;
    }
	
    // A function to print the entire network. Used for debugging
    void printNetwork() {
        cout << "Bayesian Network Structure: " << endl;
        int index = 0;
        list<Graph_Node>::iterator listIt;

        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            cout << "Node " << index << ": " << endl;
            listIt->printNode();
            cout << "------------------------------------------------------------------------------------------" << endl;
            index++;
        }
    }

    void setParentNumValues() {
        for (auto listIt = Pres_Graph.begin(); listIt!=Pres_Graph.end(); listIt++) {
            vector<string> parents = listIt->get_Parents();
            vector<int> parent_num_value;
            for (auto parent : parents) {
                parent_num_value.push_back(search_node(parent)->get_nvalues());
            }
            parent_num_values.push_back(parent_num_value);
        }
    }

    // The following function will initiate the cpts to some values (Let us say we assume uniform distribution for initialization)
    void initializeCPTs() {
        list<Graph_Node>::iterator listIt;
        int index = 0;

        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            int num_values = listIt->get_nvalues();
            int num_parents = listIt->get_Parents().size();

            int cpt_size = listIt->get_CPT().size();

            // Initialize CPT with uniform probabilities
            vector<double> cpt(cpt_size, 1.0/cpt_size);

            listIt->set_CPT(cpt);

            index++;
        }
    }

    void initializeRandomCPT() {
        list<Graph_Node>::iterator listIt;
        // Create a random number generator
        std::random_device rd; // Used to seed the random engine
        std::mt19937 gen(rd()); // Mersenne Twister 19937 generator
        std::uniform_real_distribution<double> dis(0.0, 1.0); // Uniform distribution between 0 and 1

        // Generate a random number between 0 and 1
        double randomValue = dis(gen);

        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            vector<double> curr_cpt = listIt->get_CPT();
            int cpt_size = curr_cpt.size();

            vector<double> new_cpt(cpt_size);

            for (int i = 0; i < cpt_size; i++) {
                if (curr_cpt[i] == -1) {
                    double newProb = dis(gen);
                    new_cpt[i] = newProb;
                }
                else {
                    new_cpt[i] = curr_cpt[i];
                }
            }

            listIt->set_CPT(new_cpt);
        }
    }

    // Initializes a list of list of indexes used to find the cpt index for a particular variable in the cpt.
    void initializeIndexes() {
        int net_size = this->netSize();
        cpt_indexes.resize(net_size);
        list<Graph_Node>::iterator listIt;
        int index = 0;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            vector<int> curr_indexes;
            vector<string> parents = listIt->get_Parents();
            // vector<int> parent_num_values;

            // for (int i = 0; i < parents.size(); i++) {
            //     parent_num_values.push_back(search_node(parents[i])->get_nvalues());
            // }

            int tmp = 1;

            for (int i = 0; i < parents.size(); i++) {
                tmp = tmp*(parent_num_values[index][i]);
            }

            curr_indexes.push_back(tmp);

            for (int i = 0; i < parents.size(); i++) {
                tmp = tmp/(parent_num_values[index][i]);
                curr_indexes.push_back(tmp);
            }

            cpt_indexes[index] = curr_indexes;
            index++;
        }
    }

    // A function to print the entire CPT
    void printCPTs() {
        list<Graph_Node>::iterator listIt;

        for (listIt == Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            listIt->printCPT();
        }
    }

    // A function to calculate the probability of a variable given its parents.
    double computeCondProbGivenParents(list<Graph_Node>::iterator listIt, int list_index,int node_index, vector<int>& parent_values) {
        list<Graph_Node>::iterator listIt2 = get_nth_node(node_index);
        int n = listIt->get_Parents().size();
        // int list_index = get_index(listIt->get_name());
        vector<int> cpt_index = cpt_indexes[list_index];
        int index = 0;

        for (int i = 0; i < n; i++) {
            index = index + (parent_values[i])*(cpt_index[i+1]);
        }
        index = index +  node_index*cpt_index[0];
        return listIt->get_CPT()[index];
    }

    // A function to compute probability
    double computeCondProbability(vector<int>& records) {
        double result = 1.0;
        vector<int> parent_values;

        list<Graph_Node>::iterator listIt;
        int list_index = 0;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            parent_values.clear();
            int i = get_index(listIt->get_name());

            vector<string> parents = listIt->get_Parents();
            for (string &parent : parents) {
                parent_values.push_back(records[get_index(parent)]);
            }

            int k = records[i];
            result *= computeCondProbGivenParents(listIt, list_index,k, parent_values);
            list_index++;
        }
        return result;
    }

    void initializeExpectationValues(vector<vector<double>>& expectation_values, vector<vector<int>>& patient_records, vector<bool>& has_missing_variable) {
        int num_patients = patient_records.size();
        int num_variables = patient_records[0].size();

        vector<double> curr_expectation_values;
        
        for (int i = 0; i < num_patients; i++) {
            curr_expectation_values.clear();

            if (!has_missing_variable[i]) {
                curr_expectation_values.push_back(-1);
                continue;
            }

             for (int j = 0; j < num_variables; j++) {
                int curr_value = patient_records[i][j], num_values;

                if (curr_value != -1) continue;

                num_values = get_nth_node(j)->get_nvalues();

                for (int k = 0; k < num_values; k++) {
                    curr_expectation_values.push_back(0);
                }
             }
             expectation_values.push_back(curr_expectation_values);
        }
    }

    void expectationStep(vector<vector<double>>& expectation_values, vector<vector<int>>& data_records, const int patient_number, const int missing_variable) {
        list<Graph_Node>::iterator listIt = get_nth_node(missing_variable);
        int num_values = listIt->get_nvalues();
        vector<double> new_expectations;
        double sum = 0.0;
        // Calculate the denominator
        for (int i = 0; i < num_values; i++) {
            data_records[patient_number][missing_variable] = i;
            double prob = computeCondProbability(data_records[patient_number]);
            new_expectations.push_back(prob);
            sum += prob;
            data_records[patient_number][missing_variable] = -1;
        }

        // Normalize the probabilities by dividing by the sum
        for (int i = 0; i < num_values; i++) {
            new_expectations[i] /= sum;
        }

        expectation_values[patient_number] = new_expectations;
    }

    void maximizationStep(vector<vector<double>>& expectation_values, vector<vector<int>>& data_records, double& max_change_threshold) {
        // Loop through each node in the bayesian network
        list<Graph_Node>::iterator listIt;

        vector<double> curr_cpt;

        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            // Get the cpt for the current node
            int curr_index = this->get_index(listIt->get_name());
            curr_cpt = listIt->get_CPT();

            // Find the number of values in the CPT and the number of values the current variable can take
            int cpt_size = curr_cpt.size();
            int num_values = listIt->get_nvalues();

            // Create vectors to store cpt_indexes, sum of missing values, and sum of non-missing values
            vector<int> parent_cpt_indexes;
            double total_multipliers = cpt_indexes[curr_index][0];

            vector<double> numerator(cpt_size, 0);
            vector<double> denominator(cpt_size, 0);

            // Determine the number of patients and variables in the data
            int num_patients = expectation_values.size();
            int num_variables = expectation_values[0].size();


            // loop through each record in the data
            for (int i = 0; i < num_patients; i++) {
                // Vector to store the values of the current record
                int missing_variable_index_1 = -1, missing_variable_index_2 = -1, missing_parent_index_1 = -1, missing_parent_index_2 = -1;
                vector<int> curr_patient;

                // variable to store the current node's value
                int curr_node_value = -1;

                // Get the value of the current node in the current record
                int curr_value = data_records[i][curr_index];
                // check if the current value is missing
                if (curr_value == -1) {
                    missing_variable_index_1 = 0;
                    missing_parent_index_1 = curr_index;
                }

                // Store the current node's value
                curr_node_value = curr_value;

                // Loop through parent nodes of the current node
                vector<string> parents = listIt->get_Parents();
                for (int j = 0; j < parents.size(); j++) {
                    curr_value = data_records[i][get_index(parents[j])];

                    if (curr_value == -1) {
                        missing_parent_index_1 = missing_parent_index_2 = get_index(parents[j]);
                        missing_variable_index_1 = j+1;
                        missing_variable_index_2 = j;
                    }

                    curr_patient.push_back(curr_value);
                }

                // Calculate sums for both missing and non-missing values
                computeNumerator(listIt, expectation_values, curr_cpt, numerator, i, num_values, curr_patient, missing_variable_index_1, missing_parent_index_1, curr_node_value);
                computeDenominator(listIt, expectation_values, curr_cpt, denominator, i, num_values, curr_patient, missing_variable_index_2, missing_parent_index_2, curr_node_value);

            }

            updateCPTValues(curr_cpt, numerator, denominator, num_values, max_change_threshold);
            listIt->set_CPT(curr_cpt);
        }
    }

    // function to calculate the sum of missing values for the current node's cpt
    void computeNumerator(list<Graph_Node>::iterator &listIt, vector<vector<double>>& expectations,vector<double>& curr_cpt, vector<double>& numerator, int i, int num_values, vector<int>& curr_patient, int missing_variable_index_1, int missing_parent_index_1, int curr_node_value) {
        if (missing_variable_index_1 == -1) { // no missing nodes for the current node
            int node_index = get_index(listIt->get_name());
            int num_parents = listIt->get_Parents().size();
            int index = 0;
            for (int parent_ind = 0; parent_ind < num_parents; parent_ind++) {
                index = index + (curr_patient[parent_ind])*(cpt_indexes[node_index][parent_ind + 1]);
            }

            index = index + curr_node_value * cpt_indexes[node_index][0];
            numerator[index] += 1;
        }

        else {
            int node_index = get_index(listIt->get_name());
            int num_parents = listIt->get_Parents().size();

            int num_values_miss = get_nth_node(missing_parent_index_1)->get_nvalues();

            double increment = 0.0;
            for (int j = 0; j < num_values_miss; j++) {
                increment = expectations[i][j];
                int index = 0;
                if (missing_variable_index_1 == 0) {
                    curr_node_value = j;
                }
                else {
                    curr_patient[missing_variable_index_1-1] = j;
                }

                for (int k = 0; k < num_parents; k++) {
                    index = index + (curr_patient[k])*(cpt_indexes[node_index][k+1]);
                }
                index = index + curr_node_value * cpt_indexes[node_index][0];
                numerator[index] += increment;

                if (missing_variable_index_1==0) {
                    curr_node_value = -1;
                }
                else {
                    curr_patient[missing_variable_index_1 - 1] = -1;
                }
            }
        }
    }

    void computeDenominator(list<Graph_Node>::iterator &listIt, vector<vector<double>>& expectations, vector<double>& curr_cpt, vector<double>& denominator, int i, int num_values, vector<int>& curr_patient, int missing_variable_index_2, int missing_parent_index_2, int curr_node_value) {
        if (missing_variable_index_2 == -1) {
            int node_index = get_index(listIt->get_name());
            int num_parents = listIt->get_Parents().size();
            int index = 0;

            for (int parent_index = 0; parent_index < num_parents; parent_index++) {
                index = index + (curr_patient[parent_index])*cpt_indexes[node_index][parent_index+1];
            }
            for (int j = 0; j < num_values; j++) {
                int index_1 = index + j*cpt_indexes[node_index][0];
                denominator[index_1] += 1;
            }
        }
        else {
            int num_parents = listIt->get_Parents().size();
            int node_index = get_index(listIt->get_name());

            int curr_num_value = get_nth_node(missing_parent_index_2)->get_nvalues();
            double increment = 0.0;

            for (int j = 0; j < curr_num_value; j++) {
                increment = expectations[i][j];
                int index = 0;

                curr_patient[missing_variable_index_2] = j;

                for (int p = 0; p < num_parents; p++) {
                    index = index + (curr_patient[p])*(cpt_indexes[node_index][p+1]);
                }

                for (int k = 0; k < num_values; k++) {
                    int index1 = index + k*cpt_indexes[node_index][0];
                    denominator[index1] += 1.095*increment;
                }

                curr_patient[missing_variable_index_2] = -1;
            }
        }
    }

    void updateCPTValues(vector<double>& curr_cpt, vector<double>& numerator, vector<double>& denominator, int num_values, double &max_change_threshold) {
        for (int i = 0; i < curr_cpt.size(); i++) {
            double prob = (numerator[i] +0.0001)/(denominator[i] + 0.001*num_values);
            max_change_threshold = max(max_change_threshold, fabs(curr_cpt[i] - prob));
            curr_cpt[i] = prob;

            if (curr_cpt[i] < 0.0001) { curr_cpt[i] = 0.0001;}
            if (curr_cpt[i] > 0.9999) {curr_cpt[i] = 0.9999;}
        }
    }

    void saveSolvedBayesNet(const string& outputFile) {
        ofstream solved_alarm(outputFile);

        string line, temp;

        if (!solved_alarm.is_open()) {
            cerr << "Error: Cannot save the solved network into 'solved_alarm.bif' file. Try Again!!!!" << endl;
            exit(2);
        }
        else {
            // solved_alarm << fixed;
            // solved_alarm.precision(4);

            ifstream unsolved_alarm("alarm.bif");
            if ( unsolved_alarm.is_open() ) {
                while (!unsolved_alarm.eof()) {
                    stringstream ss;
                    getline(unsolved_alarm, line);

                    ss.str(line);
                    ss >> temp;

                    if (temp.compare("probability") == 0) {
                        ss >> temp;
                        ss >> temp;

                        list<Graph_Node>::iterator listIt;
                        list<Graph_Node>::iterator listIt1;

     				    listIt=this->search_node(temp);
                        int index=this->get_index(temp);

                        solved_alarm << line << endl;

                        getline(unsolved_alarm, line);
                        solved_alarm << "\ttable";

                        vector<double> cpt = listIt->get_CPT();
                        
                        for (int i = 0; i < cpt.size(); i++) {
                            if (!unsolved_alarm.eof()) { solved_alarm << " " << fixed << setprecision(4) << cpt[i]; }
                        }
                        if (!unsolved_alarm.eof()) { solved_alarm << " ;" << endl; }
                    }

                    else {
                        if (!unsolved_alarm.eof()) {solved_alarm << line << endl; }
                    }
                }
                unsolved_alarm.close();
            }
            else {
                cerr << "Error: Unable to open alarm.bif file." << endl;
                exit(1);
            }
            solved_alarm.close();
        }
    }

};

network read_network()
{
	network Alarm;
	string line;
	int find=0;
  	ifstream myfile("alarm.bif"); 
  	string temp;
  	string name;
  	vector<string> values;
  	
    if (myfile.is_open())
    {
    	while (! myfile.eof() )
    	{
    		stringstream ss;
      		getline (myfile,line);
      		
      		
      		ss.str(line);
     		ss>>temp;
     		
     		
     		if(temp.compare("variable")==0)
     		{
                    
     				ss>>name;
     				getline (myfile,line);
                   
     				stringstream ss2;
     				ss2.str(line);
     				for(int i=0;i<4;i++)
     				{
     					
     					ss2>>temp;
     					
     					
     				}
     				values.clear();
     				while(temp.compare("};")!=0)
     				{
     					values.push_back(temp);
     					
     					ss2>>temp;
    				}
     				Graph_Node new_node(name,values.size(),values);
     				int pos=Alarm.addNode(new_node);

     				
     		}
     		else if(temp.compare("probability")==0)
     		{
                    
     				ss>>temp;
     				ss>>temp;
     				
                    list<Graph_Node>::iterator listIt;
                    list<Graph_Node>::iterator listIt1;
     				listIt=Alarm.search_node(temp);
                    int index=Alarm.get_index(temp);
                    ss>>temp;
                    values.clear();
     				while(temp.compare(")")!=0)
     				{
                        listIt1=Alarm.search_node(temp);
                        listIt1->add_child(index);
     					values.push_back(temp);
     					
     					ss>>temp;

    				}
                    listIt->set_Parents(values);
    				getline (myfile,line);
     				stringstream ss2;
                    
     				ss2.str(line);
     				ss2>> temp;
                    
     				ss2>> temp;
                    
     				vector<double> curr_CPT;
                    string::size_type sz;
     				while(temp.compare(";")!=0)
     				{
     					curr_CPT.push_back(atof(temp.c_str()));
     					ss2>>temp;
    				}
                    listIt->set_CPT(curr_CPT);
     		}
            else
            {}		
    	}
    	
    	if(find==1)
    	myfile.close();
  	}
  	
  	return Alarm;
}

// The following function will read the "records.dat" file and store the data in two vectors.
// Will be used during the execution.
void readRecords(network& alarm, const string& datafilename, vector<vector<int>>& data_records, vector<pair<int, int>>& missing_values, vector<bool>& has_missing) {
    ifstream record_file(datafilename);

    string curr_line;
    int patient_number = 0;

    if (record_file.is_open()) {
        while (!record_file.eof()) {
            stringstream ss;
            string temp;
            bool has_miss = false;

            getline(record_file, curr_line);

            ss.str(curr_line);

            vector<int> curr_record;
            int variable_number = 0;

            while (variable_number < alarm.netSize()) {
                ss >> temp;

                if (temp == "\"?\"") {
                    missing_values.push_back({patient_number, variable_number});
                    curr_record.push_back(-1);
                    has_miss = true;
                }
                else curr_record.push_back(alarm.get_nth_node(variable_number)->valueIndex(temp));
                variable_number++;
            }

            data_records.push_back(curr_record);
            has_missing.push_back(has_miss);
            patient_number++;
        }
        record_file.close();
    }
}

int main()
{
    chrono::milliseconds timeout(105000);
    auto startTime = high_resolution_clock::now();
	network Alarm;
	Alarm=read_network();
    Alarm.setParentNumValues();

    double convergenceThreshold = 1e-5;

    vector<vector<int>> data_records;
    vector<pair<int, int>> missing_values;

    vector<bool> has_missing_variable;

    readRecords(Alarm, "records.dat", data_records, missing_values, has_missing_variable);

    // Alarm.initializeCPTs();
    Alarm.initializeRandomCPT();
    Alarm.initializeIndexes();

    vector<vector<double>> expectation_values;

    Alarm.initializeExpectationValues(expectation_values, data_records, has_missing_variable);

    auto initializetime = high_resolution_clock::now();
    auto initduration = duration_cast<milliseconds>(initializetime - startTime);
    cout << "initialization time: " << initduration.count() << endl;

    int count = 0;
    // cout << "starting" << endl;

    bool breakloop = false;

    double max_change_threshold = 5.0;
    auto loopstart = high_resolution_clock::now();
    while (max_change_threshold > convergenceThreshold) {
        cout << count << endl;
        max_change_threshold = 0.0;
        count++;
        for (int i = 0; i < missing_values.size(); i++) {
            auto time = high_resolution_clock::now();
            auto curr_duration = duration_cast<milliseconds>(time - startTime);
            if (curr_duration > timeout) {
                // cout << "breaking " << endl;
                breakloop = true;
                break;
            }
            int patient_number, missing_variable;
            patient_number = missing_values[i].first;
            missing_variable= missing_values[i].second;
            Alarm.expectationStep(expectation_values, data_records, patient_number, missing_variable);
        }

        if (breakloop) {
            break;
        }

        Alarm.maximizationStep(expectation_values, data_records, max_change_threshold);
        cout << max_change_threshold << endl;
    }
    auto loopend = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(loopend - loopstart);
    // cout << "Time taken for EM algorithm: " << duration.count() << "seconds" << endl;

    // Alarm.printNetwork();

    Alarm.saveSolvedBayesNet("solved_alarm.bif");

    auto final_time = high_resolution_clock::now();

    auto finalduration = duration_cast<seconds>(final_time - startTime);
    cout << "Total time taken: " << finalduration.count() << "seconds" << endl;

    const char* cmd = "./format_checker";
    int exit = system(cmd);

	// cout<<"Perfect! Hurrah! \n";
}