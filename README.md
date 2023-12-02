# A4 Medical Diagnosis

## Goal
The goal of this assignment is to gain experience with learning Bayesian Networks and understanding their value in real-world applications.

## Scenario
In this assignment, you'll be working on a medical diagnosis scenario. Medical researchers have created a Bayesian network to model the relationships between diseases and observed symptoms. Your task as a computer scientist is to learn the parameters for this network based on health records. Some records have missing values, and your goal is to compute these parameters, making the network suitable for medical diagnosis.

## Problem Statement
You are provided with a Bayesian Network structure in the .bif (Bayesian Interchange Format) format, representing the relationships between diseases and symptoms. Probability values in the network are missing and represented as -1. You also have a dataset file with patient records, where each record includes observed variable values and may contain missing values marked as '?'.

## Input Format
- The input consists of the Bayesian network structure in .bif format, where most probability values are -1.
- The dataset file (e.g., records.dat) includes patient records with known and unknown variable values.

## Output Format
Your task is to learn the missing probability values and replace them with computed values, accurate up to four decimal places. Your output should be a complete .bif file representing the Bayesian network with learned parameters.

## What's Provided
- Dataset file (records.dat) with patient records.
- Starter code for parsing the .bif file and data structure representation.
- A format checker to validate your output.

## What to Submit
1. Compressed code in a .zip file named according to the specified format.
2. Writeup.txt with your choice of programming language (C++, Java, or Python), and names of collaborators if any.

## Code Verification
- Ensure that your code follows the input/output specifications provided.
- Test your code on the provided dataset (records.dat) to verify accurate parameter learning.

## Evaluation Criteria
1. Accuracy of learned parameter values.
2. Extra credit may be awarded for standout performance.

## Allowed and Not Allowed
- Work individually or in a group of two.
- Choose C++, Java, or Python.
- Do not use built-in libraries for Bayes nets or expectation maximization.
- Do not discuss the assignment with anyone outside the class, and follow academic integrity guidelines.

## Submission Reminder
- Submit a .zip file named correctly and adhering to the specified format.
- Make sure your code compiles and runs on the required machine.
- Provide a run.sh script that takes alarm.bif and a data file as inputs and generates a solved_alarm.bif output.

## Code Verification
- Your code should be able to finish learning in 2 minutes on a dataset with around 10,000 patient records.

### Diagnosis
1. Hypovolemia
2. Left Ventricular Failure
3. Anaphylaxis
4. Insufficient Analgesia
5. Pulmonary Embolus
6. Intubation
7. Kinked Tube
8. Disconnection

### Observed Variables (Symptoms):
1. CVP (Central Venous Pressure)
2. PCWP (Pulmonary Capillary Wedge Pressure)
3. History
4. TPR (Total Peripheral Resistance)
5. Blood Pressure
6. CO (Cardiac Output)
7. HR BP (Heart Rate from Blood Pressure)
8. HR EKG (Heart Rate from Electrocardiogram)
9. HR SAT (Heart Rate from saturation)
10. SaO2 (Oxygen Saturation)
11. PAP (Pulmonary Arterial Pressure)
12. MV (Mechanical Ventilation)
13. Min Vol (Minimum Ventilation)
14. Exp CO2 (Expired CO2)
15. FiO2 (Fraction of Inspired Oxygen)
16. Pres (Pressors)