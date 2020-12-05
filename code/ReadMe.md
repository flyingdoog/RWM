# Complie
  g++ ./code/*.cpp -o RWM


# Run
./RWM -d 6ng


# Parameters

  * -d data directory path
  
  * -n the target network
   
  * -a alpha: default: 0.9
  
  * -l lambda: default: 0.7
  
  * -t theta: default: 0.9
  
  * -e epsilon: default: 0.01


# Input file format

 ## Edge list
 
The edge lists of networks are store in path_to_your_data_directory/networks.txt

The first line expresses the number of networks, and the second line points out whether the networks are directed or not.

For example:
```
#networks:5
#directed:0
```

Then, edge lists of networks are followed. networks starts by '--'

For exmaple
```
--0--
/*edge lists of network 0*/
--1--
/*edge lists of network 1*/
...
--k--
/*edge lists of network k*/
```

RWM supports both weighted and unweighted networks.

For weighted networks, the input format for an edge is "source_node	target_node	weight" (delimiter is "\t"). For example:

```
1	88160	1.3
22	118052	2.1
35	161555	7.5
14	244916	0.5
25	346495	3.8
```

For unweighted graph, the input format for an edge is "source_node	target_node" (delimiter is "\t"). For example:

 ```
1	88160
22	118052
35	161555
14	244916
25	346495
```

## Cross networks edges

The edge lists of cross network graphs are store in path_to_your_data_directory/inter.txt

Each network starts with 

=i,j=

and ends with 

--end--

where i and j are networks the graph connected.

For example
```
=0,1=
/*cross networks edges from  network 0 to network 1*/
--end--
=0,2=
/*cross networks edges from  network 0 to network 2*/
....
```

Both weighted and unweighted graphs are supported. Formats are similar to the above networks.


## Community file for evaluation
The community of nodes in network $i$ are store in path_to_your_data_directory/cmty$i$.txt
with each line represents a community, splitted with '\t'


# Output
This version will not write any files to the disk. The terminal will output the accurarcy results (Macro pre, Macro Recall, Macro F1).
