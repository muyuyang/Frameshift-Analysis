from ete3 import Tree

# Oryza
def generate_newick():
    in1 = Tree()
    in2 = in1.add_child()
    oc = in2.add_child(name='Oryza coarctata')
    in3 = in2.add_child()
    in4 = in3.add_child()
    in5 = in3.add_child()
    in6 = in4.add_child()
    in7 = in4.add_child()
    in8 = in6.add_child()
    in9 = in6.add_child()
    in10 = in8.add_child()
    ogl = in8.add_child(name='Oryza glaberrima')
    os = in10.add_child(name='Oryza sativa')
    in11 = in10.add_child()
    on = in11.add_child(name='Oryza nivara')
    orr = in11.add_child(name='Oryza rufipogon')

    op = in9.add_child(name='Oryza punctata')
    om = in9.add_child(name='Oryza minuta')

    oo = in7.add_child(name='Oryza officinalis')
    oa = in7.add_child(name='Oryza alta')

    oau = in5.add_child(name='Oryza australiensis')
    
    in13 = in1.add_child()
    ob = in13.add_child(name='Oryza brachyantha')
    orri = in13.add_child(name='Oryza ridleyi')

    print(in1.write())




generate_newick()