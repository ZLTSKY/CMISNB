1. Download the initial network from STRING；
2. Match with cancer data and delete the edge without data from the network;
3. Remove the repeated edges to get the initial background network;
4. background.py is used to obtain a background network that can be used for CMISNB calculation;
5. Input reference sample data, cancer sample data and background network into knncmi_network.py to calculate CMI.