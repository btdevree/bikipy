#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create an example network!
"""
import networkx as nx
import bikipy.bikicore.model as bkcm
import bikipy.bikicore.components as bkcc

# Make a new model
m1 = bkcm.Model(1, 'dDAT Model', None)

## Make the component pieces
#d1 = bkcc.Drug()
#d1.name = 'Dopamine'
#d1.symbol = 'DA'

d2 = bkcc.Drug()
d2.name = 'Sodium ion 1'
d2.symbol = 'Na1'

d3 = bkcc.Drug()
d3.name = 'Sodium ion 2'
d3.symbol = 'Na2'

d4 = bkcc.Drug()
d4.name = 'Chloride ion'
d4.symbol = 'Cl'

d5 = bkcc.Drug()
d5.name = '[3H]-Beta-CFT'
d5.symbol = '[3H]B'

p1 = bkcc.Protein()
p1.name = 'Drosophila dopamine transporter'
p1.symbol = 'DAT'
p1.conformation_names = ['Outward Open', 'Outward Closed', 'Inward Closed', 'Inward Open']
p1.conformation_symbols = ['OO', 'OC', 'IC', 'IO']

#m1.drug_list.extend([d1, , d5])
m1.drug_list = [d2, d3, d4, d5]
m1.protein_list = [p1]

## Make the rules
#r1 = bkcc.Rule(m1)
#r1.rule_subject = [d1]
#r1.subject_conf = [None]
#r1.rule = ' reversibly associates with '
#r1.rule_object = [p1]
#r1.object_conf = [[0]]

r2 = bkcc.Rule(m1)
r2.rule_subject = [d2]
r2.subject_conf = [None]
r2.rule = ' reversibly associates with '
r2.rule_object = [p1]
r2.object_conf = [[0]]

r3 = bkcc.Rule(m1)
r3.rule_subject = [d3]
r3.subject_conf = [None]
r3.rule = ' reversibly associates with '
r3.rule_object = [p1]
r3.object_conf = [[0]]

r4 = bkcc.Rule(m1)
r4.rule_subject = [d4]
r4.subject_conf = [None]
r4.rule = ' reversibly associates with '
r4.rule_object = [p1]
r4.object_conf = [[0]]

#r5 = bkcc.Rule(m1)
#r5.rule_subject = [d1]
#r5.subject_conf = [None]
#r5.rule = ' reversibly associates with '
#r5.rule_object = [p1]
#r5.object_conf = [[3]]

r6 = bkcc.Rule(m1)
r6.rule_subject = [d2]
r6.subject_conf = [None]
r6.rule = ' reversibly associates with '
r6.rule_object = [p1]
r6.object_conf = [[3]]

r7 = bkcc.Rule(m1)
r7.rule_subject = [d3]
r7.subject_conf = [None]
r7.rule = ' reversibly associates with '
r7.rule_object = [p1]
r7.object_conf = [[3]]

r8 = bkcc.Rule(m1)
r8.rule_subject = [d4]
r8.subject_conf = [None]
r8.rule = ' reversibly associates with '
r8.rule_object = [p1]
r8.object_conf = [[3]]

r9 = bkcc.Rule(m1)
r9.rule_subject = [d5]
r9.subject_conf = [None]
r9.rule = ' reversibly associates with '
r9.rule_object = [p1]
r9.object_conf = [[0]]

r10 = bkcc.Rule(m1)
r10.rule_subject = [d5]
r10.subject_conf = [None]
r10.rule = ' reversibly associates with '
r10.rule_object = [p1]
r10.object_conf = [[1]]

r11 = bkcc.Rule(m1)
r11.rule_subject = [p1]
r11.subject_conf = [[0]]
r11.rule = ' reversibly converts to '
r11.rule_object = [p1]
r11.object_conf = [[1]]

r12 = bkcc.Rule(m1)
r12.rule_subject = [p1]
r12.subject_conf = [[1]]
r12.rule = ' reversibly converts to '
r12.rule_object = [p1]
r12.object_conf = [[2]]

r13 = bkcc.Rule(m1)
r13.rule_subject = [p1]
r13.subject_conf = [[2]]
r13.rule = ' reversibly converts to '
r13.rule_object = [p1]
r13.object_conf = [[3]]

#m1.rule_list = [r1, r2, r3, r4, r5, r6, r7, r8]
m1.rule_list = [r2, r3, r4, r6, r7, r8, r9, r10, r11, r12, r13]
#m1.rule_list = [r5, r6, r7, r8]
# Make model graph
m1.generate_network(save_graphs=True)

# Print graph
m1._fancy_main_graph_draw('testgraph')