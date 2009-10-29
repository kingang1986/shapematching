#!/home/longbin/install/apps/Python-2.5_gcc/bin/python

import sys
import math

import getopt
import string
import fileinput

clrs  = ['ff0000', '00ff00', '0000ff']


def printheader(n):
    print '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
    print '   <graphml xmlns="http://graphml.graphdrawing.org/xmlns/graphml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/graphml http://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd">'
    print '<key for="node" id="d0" yfiles.type="nodegraphics"/>'
    print '<key attr.name="description" attr.type="string" for="node" id="d1"/>'
    print '<key for="edge" id="d2" yfiles.type="edgegraphics"/>'
    print '<key attr.name="description" attr.type="string" for="edge" id="d3"/>'
    print '<key for="graphml" id="d4" yfiles.type="resources"/>'
    print '<graph edgedefault="directed" id="G" parse.edges="1" parse.nodes="%d" parse.order="free">'%(n)
def printtail():
    print '</graph>'
    print '<data key="d4">'
    print '<y:Resources/>'
    print '</data>'
    print '</graphml>'

def printPts(X, Y, s, M):
    for i in range(len(X)):
        print '    <node id="%s%d">'%(s, i)
        print '      <data key="d0">'
        print '        <y:ShapeNode>'
        print '          <y:Geometry height="16.0" width="16.0" x="%f" y="%f"/>'%(X[i], Y[i])
        if (M[i] == 0):
            print '          <y:Fill color="#000000" transparent="false"/>'
        else:
            if (s == 'n'):
                print '          <y:Fill color="#FFCC00" transparent="false"/>'
            else:
                print '          <y:Fill color="#FFCC00" transparent="false"/>'
        print '          <y:BorderStyle color="#000000" type="line" width="1.0"/>'
        print '          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="18.701171875" modelName="internal" modelPosition="c" textColor="#000000" visible="true" width="11.0" x="2.0" y="-1.8505859375">%d</y:NodeLabel>'%(i + 1)
        print '          <y:Shape type="ellipse"/>'
        print '        </y:ShapeNode>'
        print '      </data>'
        print '      <data key="d1"/>'
        print '    </node>'

def printWeightPts(X, Y, s, M, weights):
    for i in range(len(X)):
        print '    <node id="%s%d">'%(s, i)
        print '      <data key="d0">'
        print '        <y:ShapeNode>'
        print '          <y:Geometry height="%f" width="%f" x="%f" y="%f"/>'%(weights[i], weights[i], X[i], Y[i])
        if (M[i] == 0):
            print '          <y:Fill color="#000000" transparent="false"/>'
        else:
            if (s == 'n'):
                print '          <y:Fill color="#FFCC00" transparent="false"/>'
            else:
                print '          <y:Fill color="#FFCC00" transparent="false"/>'
        print '          <y:BorderStyle color="#000000" type="line" width="1.0"/>'
        print '          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="18.701171875" modelName="internal" modelPosition="c" textColor="#000000" visible="true" width="11.0" x="2.0" y="-1.8505859375">%d</y:NodeLabel>'%(i + 1)
        print '          <y:Shape type="ellipse"/>'
        print '        </y:ShapeNode>'
        print '      </data>'
        print '      <data key="d1"/>'
        print '    </node>'

def printEdge(A, B):
    for i in range(len(A)):
        print '    <edge id="e%d" source="n%d" target="m%d">'%(i, A[i], B[i])
        print '      <data key="d2">'
        print '        <y:PolyLineEdge>'
        print '          <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0"/>'
        print '          <y:LineStyle color="#%s" type="solid" width="1.0"/>'%(clrs[ 2])
        print '          <y:Arrows source="short" target="short"/>'
        print '          <y:EdgeLabel alignment="center" distance="2.0" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="4.0" modelName="six_pos" modelPosition="tail" preferredPlacement="anywhere" ratio="0.5" textColor="#000000" visible="true" width="4.0" x="90.55342102050781" y="13.346414794921884"/>'
        print '          <y:BendStyle smoothed="false"/>'
        print '        </y:PolyLineEdge>'
        print '      </data>'
        print '      <data key="d3"/>'
        print '    </edge>'

if __name__ == '__main__':
    # load the image gived on the command line
    X = []
    X1 = []
    Y = []
    Y1 = []
    U = []
    V = []
    for line in fileinput.input(sys.argv[1]):
        dr = line.split(' ')
        X.append(string.atof(dr[0]) )
        Y.append(string.atof(dr[1]))
    for line in fileinput.input(sys.argv[2]):
        dr = line.split(' ')
        X1.append(string.atof(dr[0]) )
        Y1.append(string.atof(dr[1]))
    for line in fileinput.input(sys.argv[3]):
        dr = line.split(' ')
        U.append(string.atof(dr[0]) )
        V.append(string.atof(dr[1]))


    M = []
    for i in range(len(X)):
        M.append(1)
    printheader(len(X))
    printPts(X, Y,'n', M)
    printPts(X1, Y1,'m', M)
    printEdge(U, V)
    
    printtail()
