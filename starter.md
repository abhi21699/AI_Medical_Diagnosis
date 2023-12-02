# Bayesian Network Starter Code Documentation

This document provides documentation for the C++ code that represents a Bayesian network. The code consists of classes and functions for creating, manipulating, and analyzing Bayesian networks.

## `Graph_Node` Class
- This class represents a node in a Bayesian network.
- **Private Members**:
  - `Node_Name`: A string representing the name of the node (variable).
  - `Children`: A vector of integers that stores the indices of child nodes in the graph.
  - `Parents`: A vector of strings that stores the names of the parent nodes.
  - `nvalues`: An integer representing the number of categories (possible values) the variable can take.
  - `values`: A vector of strings representing the categories or possible values of the variable.
  - `CPT`: A vector of floats representing the conditional probability table for the node.
- **Public Member Functions**:
  - `Graph_Node(string name, int n, vector<string> vals)`: Constructor to initialize the node with a name and its possible values.
  - `string get_name()`: Returns the name of the node.
  - `vector<int> get_children()`: Returns the indices of child nodes.
  - `vector<string> get_Parents()`: Returns the names of parent nodes.
  - `vector<float> get_CPT()`: Returns the conditional probability table.
  - `int get_nvalues()`: Returns the number of possible values.
  - `vector<string> get_values()`: Returns the possible values of the variable.
  - `void set_CPT(vector<float> new_CPT)`: Sets the conditional probability table.
  - `void set_Parents(vector<string> Parent_Nodes)`: Sets the parent nodes.
  - `int add_child(int new_child_index)`: Adds a child node to the list of children.
  - `void printNode()`: Prints information about the node, including its name, values, parents, children, and CPT.

## `network` Class
- This class represents the entire Bayesian network, which is a collection of `Graph_Node` objects.
- **Private Member**:
  - `Pres_Graph`: A list of `Graph_Node` objects representing the network.
- **Public Member Functions**:
  - `int addNode(Graph_Node node)`: Adds a `Graph_Node` to the network.
  - `int netSize()`: Returns the number of nodes in the network.
  - `int get_index(string val_name)`: Gets the index of a node with a given name.
  - `list<Graph_Node>::iterator get_nth_node(int n)`: Gets the iterator of a node at the nth index.
  - `list<Graph_Node>::iterator search_node(string val_name)`: Gets the iterator of a node with a given name.
  - `void printNetwork()`: Prints information about the entire Bayesian network, including the nodes' details.
  - `vector<string> getVariableNames()`: Returns a vector of strings containing the names of all the variables in the network.

## `read_network()` Function
- This function reads a Bayesian network from a file named "alarm.bif" and populates the `network` object.
- It uses the information in the file to create `Graph_Node` objects and establishes parent-child relationships between nodes based on conditional probability tables.
- Returns a `network` object representing the Bayesian network.

## `main()` Function
- The main function creates a `network` object named "Alarm" by calling the `read_network()` function.
- It can be used to interact with the Bayesian network, analyze its structure, and perform various operations.

**Note**: The code is designed to read Bayesian network information from a specific file format ("alarm.bif") and construct a Bayesian network representation in memory. It provides basic functionalities to access and manipulate the network structure.
