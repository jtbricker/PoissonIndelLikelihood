#!usr/bin/python

from collections import namedtuple
from sys import stdout

Node = namedtuple('Node', 'name, dist, arr, left, right')

leaf_probs = {'A':[1,0,0,0],'T':[0,1,0,0],'C':[0,0,1,0], 'G':[0,0,0,1]}
def_arr = [-1,-1,-1,-1]

tree = Node('O',None,def_arr,
            Node('N',20,def_arr,
                 Node('M',30,def_arr,
                      Node('A',50,leaf_probs['A'],None,None),
                      Node('I',55, def_arr,
                           Node('B',60,leaf_probs['A'],None,None),
                           Node('C',65,leaf_probs['A'],None,None)
                           )
                      ),
                 Node('J',35,def_arr,
                      Node('D',45,leaf_probs['T'],None,None),
                      Node('E',40,leaf_probs['C'],None,None)
                      )
                 ),
            Node('L',70,def_arr,
                 Node('K',80,def_arr,
                      Node('F',85,leaf_probs['C'],None,None),
                      Node('G',90,leaf_probs['G'],None,None)
                      ),
                 Node('H',75,leaf_probs['G'],None,None)
                 )
            )

tree2 = Node('5',None,def_arr,
                Node('4',15,def_arr,
                     Node('1',5,leaf_probs['A'],None,None),
                     Node('2',10,leaf_probs['T'],None,None)
                     ),
                Node('3',20,leaf_probs['C'],None,None)
            )

def printwithspace(i):
    stdout.write("%s " %i)


        
        return likelihood(node)

def likelihood(node):
    if node.left == None and node.right == None:
        return node
    else:
        for i in range(4):
            left = 0
            right = 0
            for j in range(4):
                left += substitution_prob(i,j,node.left.dist)*node.left.arr[j]
                right += substitution_prob(i,j,node.right.dist)*node.right.arr[j]
            node.arr[i] = left*right
        return node

def substitution_prob(i,j,t):
    if i == j: #no net substitution
        return 1
    elif i+j == 3: #transition (assuming A=0, T=1, C=2, G=3
        return -t
    else: #transversion
        return t

preorder(tree2)
print tree2

