#import "@preview/touying:0.5.3": *
#import "@preview/clean-math-presentation:0.1.0": *

#show: clean-math-presentation-theme.with(
  config-info(
    title: [AlfaPang: Alignment Free for pangenome construction],
    short-title: [AlfaPang: Bioinformatics seminar \@ Unipi],
    authors: (
      (name: "Luigi di Micco", affiliation-id: "1, 2"),
      // (name: "Second Author", affiliation-id: 2),
      // (name: "Third Author", affiliation-id: 1)
    ),
    author: "Luigi di Micco",
    affiliations: (
      (id: 1, name: "CS Department, University of Pisa"),
      (id: 2, name: "I am not the author of the work though"),
    ),
    date: datetime(year: 2024, month: 11, day: 20),
  ),
  config-common(
    slide-level: 2,
    //handout: true,
    //show-notes-on-second-screen: right,
  ),
  progress-bar: true,
)

#title-slide(
  logo1: image("images/University_of_Pisa.svg", height: 4.5em),
)

// == Outline <touying:hidden>

// #components.adaptive-columns(outline(title: none))

= Introduction

#slide(title: "Introduction")[
  Pangenome (or variation) graphs:
  
  - serve as models for joint representation of populations of genomes
  
  - useful in analyzing sequence evolution and variation
  
  - helpful in reducing the _"reference bias"_ in the analysis of experimental data

#pause
  However, the success of the pangenome-based approaches depends on the existence of *efficient construction methods*, applicable to large collections of genomes 
]

#slide(title: "Pangenome construction - existing approaches")[
  
  
  Most of the pangenome building algorithms adapt the approaches used in the whole genome alignment tools.
  
  - Early versions of the `VG` toolkit constructed pangenome graphs iteratively, i.e. aligning consecutive sequences to the current graph
  - The current version aligns all genomes to a reference genome (see `Minigraph-Cactus`)
  
In both approaches the outcome depends on an arbitrary choice: either genome order or reference.

Alternatives to avoid such biases:
- `seqwish`: construction does not scale linearly with the number of genomes, and the final graph needs refinements
$=>$ *pggb* pipeline: `wfmash` + `seqwish` + (`smoothxg`+`gfaffix`)
]

#slide(title: "de Bruijn graphs ")[
  de Bruijn graphs are one of the most important alternative graph pangenome models to the variation graphs
    
  - their structure is strictly determined by the parameter $k$
    - avoids order and reference biases
  
  - the construction is conceptually simple
  
  - optimized construction algorithms (see `TwoPaCo` or `bifrost`): orders of magnitude faster than alignment-based ones
  #pause
  They pose challenges for downstream analysis:
  - annotation
  - visualization
  - information extraction
]

= AlfaPang: the algorithm

#focus-slide[
  Thus, AlfaPang!
]
#slide(title: "AlfaPang: what's that? ")[
  
  `AlfaPang`:
  
  - algorithm building a _variation graph_ efficiently
  
  #pause
  
  Spoiler:
  
  *pggb* pipeline: 
  #alternatives[#strike[#text(red)[`wfmash` + `seqwish`]]][*`AlfaPang`*]+ (`smoothxg`+`gfaffix`)
  
  #pause
  
  $=>$ 
  
  - significant efficiency improvements
  - similar graphs' properties  
]
#slide(title: "AlfaPang: a lookahead ")[
  
  - *Inputs*: 
    - a collection of sequences $S$ 
    - a positive natural number $k$

    
- *Algorithm steps*:
  
  + Build the generic representation graph $G=<V,E>$
  
  + Build a weighted bipartite graph with parts $V$ and $B$
    - $B$ is a set of vertices labeled by canonical $k$-mers of $S$
    - each edge $e$ is assigned a value $C(e) in {-k, ..., -1,1, ..., k}$
      - $C(<<i,j>, b>)=c$ means that: 
        + the position in the sequence $<i, j>$  can be extended to a $k$-mer represented by $b$
        + $c$ indicates the position of $S_i [j]$ in the $k$-mer 
  
        
  + Build the quotient graph $G'$:
    - Traverse the bipartite graph in a BFS fashion starting from a vertex $v$
      - There are constraints to adhere to during the visit
    - Every vertex belonging to $V$, visited during one such run, establishes an equivalence class  
      // - $C(<<i,j>, b>)=c <=>$
      // $ S_i [j-c+1 .. j+k-c]=l(b)$ for $c > 0$ 
      
      // $ S_i [j-c-k .. j-c+1]=l(b)^(-1) "for" c < 0$ 
    
]
= Background concepts
#focus-slide[
  A step back
]
#slide(title: "Directed variation graphs ")[
  #linebreak()
  
  #definition(//title: "Important stuff"
)[
    A _directed variation graph_ is a tuple $G = <V, E, l>$, where:
    - $V$ is a set of vertices,
    - $E subset.eq V^2$ is a set of _directed_ edges,
    - $l: V-> Sigma^(\+)$ is a function labeling vertices with non-empty strings over the DNA alphabet $Sigma = {A, C, G, T}$
  ]
]
#slide(title: "Quotient graph ")[ 
  We assume that:
  - $G = <V, E, l>$ is a variation graph
  - $~$ is an equivalence relation on the set of G nodes, s.t.
    - $v ~ v' => l(v) = l(v')$ for all $v,v' in V$
  
  #definition(//title: "Important stuff"
)[
    The _quotient graph_ of $G "by" ~$ is defined as $G'=<V', E', l'>$, where:
    - $V'= V \/ ~$,
    - $E' = {<[v]_~, [v']_~ > | <v, v'> in E} $ is a set of _directed_ edges,
    - $l'([v]_~) = l(v)$.
  ]
  When $|l(v)|=1$ for every vertex $v$, the graph is called _singular_.
]
#focus-slide[So? 

What do we do with these?
#pause


We want to encode the sequences through them
]

#slide(title: "How to represent collections of sequences?")[
  Given: 
  - a set of sequences $cal(S)={S_1, dots, S_n}$
  - a _singular directed variation graph_ $G = <V, E, l>$
  - $pi: S -> cal(P)(G)$ 
  We say that $<G, pi>$ *represents* S iff:
  
  - $hat(l)(pi(S_i))=S_i,  forall i in {1, dots, n}$
  - every vertex in $G$ occurs in some path $pi(S_i)$
  - every edge in $G$ joins two consecutive vertices in some path $pi(S_i)$
]
#slide(title: "The generic representation of a collection of strings")[
  We define
  -  $"Pos"(cal(S))={<i,j> | 1 lt.eq i lt.eq n and 1 lt.eq j lt.eq |S_i|}$
  
  Let's call the *generic representation* of $cal(S)$ the representation $<G_0, pi_0>$ where $G = <V_0, E_0, l_0>$, and:
  
  - $V_0 = "Pos"(cal(S))$
  - $E_0={ <<i,j>, <i,j+1>> | 1 lt.eq i lt.eq n and 1 lt.eq j lt |S_i|}$
  - $l_0(<i,j>)=S_i [j]$
  - $pi_0(S_i)=<<i,1>, dots,<i,|S_i|>>$
  #align(horizon)[
    #figure(
     image("images/generic_vg.png"),
    caption: [Generic representation variation graph for the strings on the left. ])
  ]
]
#slide(title: "Sequences representation: a few more theoretical guarantees")[
  #linebreak()
  #lemma()[For every singular representation $<G, pi>$ of $cal(S)={S_1, dots, S_n}$ 
  there exist equivalence relation $~_(<G,pi>)$ on $"Pos"(cal(S))$ such that:
  + $G$ is isomorphic to a quotient graph $G_0$ by $~_(<G,pi>)$
  + $pi(S_i)=<[<i,1>]_(~_(<G,pi>)), dots, [<i,|S_i|>]_(~_(<G,pi>))> forall i in {1, dots, n}$
  ]
  - $~_(<G,pi>)$ is defined s.t.: $<i,j> ~_(<G,pi>) <i',j'> <=> pi(S_i)[j] = pi(S_(i'))[j']$
]

#slide(title: "Sequences representation: k-mer")[
  #linebreak()
  Let $S_i, S_(i') in cal(S)$, assume that they have a common $k$-mer $S_i [p dots p+k-1] = S_(i')[p' dots p'+k-1]$
  #linebreak()
  
  A few useful definitions:
  - $pi$ _reflects_ this common $k$-mer $<=>pi(S_i) [p dots p+k-1] = pi(S_(i'))[p' dots p'+k-1] $;
  - $<G,pi>$ _represents_ $cal(S) " "k"-completely" <=> $ all common $k$-mers are reflected by $pi$;
  Furthermore:

  We say the pair of $pi$-occurrences $<i,j>, <i',j'>$ of a vertex $v$ is:
  - *directly $k$-extendable* $<=>pi(S_i) [j-m dots j+m'] = pi(S_(i'))[j'-m dots j'+m'] "for some" m,m' gt.eq 0 "s.t." m+m'gt.eq k-1 $;
  - *$k$-extendable* if there's a sequence of occurrences of $v$ that starts from $<i,j>$, ends at $<i',j'>$ and each two consecutive occurrences in that sequence are directly $k$-extendable.
  We say that $G,pi>$ _represents_ $cal(S) k"-faithfully"$ if every pair of occurrences of a vertex is $k$-extendable.
  Eventually, we have:
  #theorem()[Let $cal(S)={S_1, dots, S_n}$, then the k-complete and k-faithful representation of S exists and is unique up to isomorphism.

  ]
  - Proof key fact: definition of the binary relation $~_0 $ indicating the pairs of positions in $cal(S)$ that should be merged in a representation reflecting commonn $k$-mers:
  
  $<i,j> ~_0 <i',j'> <=> exists_(0 lt.eq m < k) S_i [j-m dots j+k-1-m] = S_(i')[j'-m dots j'+k-1-m]$ 
  
]
#focus-slide()[
  People waiting for the algorithm details:
  #image("images/meme.png")
]
= AlfaPang algorithm

#slide(title: "(back to) AlfaPang")[
    - *Inputs*: a collection of sequences $S$ and a positive natural number $k$
    - *Outputs*: the quotient graph $G'$
    
- *Algorithm steps*:
  
  + Build the generic representation graph $G=<V,E>$ #emoji.checkmark.box
  
  + Build a weighted bipartite graph with parts $V$ and $B$
    - $B$ is a set of vertices labeled by canonical $k$-mers of $S$ #emoji.checkmark.box
    - each edge $e$ is assigned a value $C(e) in {-k, ..., -1,1, ..., k}$ #emoji.magnify.r
      - $C(<<i,j>, b>)=c$ means that: 
        + the position in the sequence $<i, j>$  can be extended to a $k$-mer represented by $b$
        + $c$ indicates the position of $S_i [j]$ in the $k$-mer 
]

#slide(title: "AlfaPang: edges weighting")[
  
  - $C(<<i,j>, b>)=c <=>$
      - $ S_i [j-c+1 .. j+k-c]=l(b)$ for $c > 0$
      - $ S_i [j-c-k .. j-c+1]=l(b)^(-1) "for" c < 0$ 
  
  - Embeds the concept of *bidirected variation graphs*
    - They naturally represent the double-stranded structure of DNA
    
  
    #pause
  
  3. Build the quotient graph $G'$:
    - Traverse the bipartite graph in a BFS fashion starting from a vertex $v$
      - There are _constraints to adhere_ to during the visit #emoji.magnify.r
    - Every vertex belonging to $V$, visited during one such run, establishes an equivalence class  
  
]
#slide(title: "AlfaPang: the algorithm in practice")[
  In the algorithm implementation:
  - the graph is not built explicitly
  
  But, the following data structures are used:
  #image("images/alfa_ds.drawio.svg")
  
]

#slide(title: "AlfaPang: the algorithm in practice - data structures")[
  #image("images/alfa_ds_index.drawio.png")

  - |`R`|$= |{"distinct canonical k-mers in S"}| $
    - *Inverted index*: fast locating of $k$-mer occurrences in `s`
  #pause
  - `F`: maps positions of `s` to vertex IDs of the output graph
]

#slide(title: "AlfaPang: finally some code!")[
  #image("images/quotient_code.png")

]
#slide(title: "AlfaPang: a run example")[
  We want to move from the bipartite graph to the quotient graph:
  
  #image("images/bipartite_quotient_graphs.png")
  Now, say we want to find the *#text(fuchsia)[pink]* equivalence class.

]
#slide(title: "AlfaPang: a run example")[
  #image("images/ds_beginning.png")
  - We have that `F` at that point is:
  `F = [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]`
  and `v=2`
  #pause 
  - To do that, let's start with the symbol "A" at position 3, and so:
    - `F[3] = 2` (v)
    - `K[3]` $->$ *2* $=>$ `R[`*2*`] = [3, 15]`
]
#slide(title: "AlfaPang: a run example")[
  - Next, we backtrack one position to `K[2]`
      - `K[2]` $->$ *1* $=>$ `R[`*1*`] = [2`#emoji.checkmark.box, *9*`]`
      - $=>$ `F[10]=2`(9+1)
      
    #pause
  
  - Finally, we backtrack one more position to `K[1]`
    - but `K[1]`$->$ #text(red)[0] $=>$ not the third symbol in any $k$-mer
  - 10 and 15 were pushed into the queue, so the procedure is repeated for both of them
    - No additional positions identified though
]
#slide(title: "AlfaPang: a run example")[
  So, we end up having: 
  
  `F = [0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0]`
  
  `s = [$, G, A, T, G, C, $, A, G, A, T, T, $, T, A, T, G, A, $]`
  #figure(
    image("images/subquotient_graph.png"),
    caption: [Quotient graph that we obtain from the current F. ])
]
= Experiments
#slide(title: "Experiments!")[
  #align(horizon)[
    #figure(
      image("images/mad-scientist.png")
    )
  
  ]
]
#slide(title: "Design of experiments: datasets and parameter setting", repeat: 3, self => [
  - Tested on two series of genome collections:
    - Escherichia coli #uncover("2-")[: 50, 100, 200 and 400 haplotypes, 250Mbp - 17Gbp total lengths]
    - Saccharomyces cerevisiae #uncover("2-")[: 16, 32, 64 and 118 haplotypes, 195Mbp - 1.44Gbp total lengths]
 #uncover("3-")[
    #align(horizon)[
     * How do we choose $k$?* 
   ]
 ]
])
#slide(title:"Some evaluation for k choice", repeat: 5, self => [
  #let (uncover, only, alternatives) = utils.methods(self)
  - Before the analysis, the complexity and repetitiveness were estimated
  - Calculation of the fractions of $k$-mers:
    - _Overrepresented_ #uncover("2-")[: occurring more frequently than the number of genomes]
    - _Rare_ #uncover("2-")[: occurring only once]
  
    #uncover("3-")[
      Too low $k =>$ large fraction of overrepresented $k$-mers #uncover("4-")[$=>$ merge of non-homologous fragments]
     
      Too large $k =>$ increase the fraction of rare $k$-mers  #uncover("4-")[$=>$ poor sequence similarity detection]
    ]
    #uncover("5-")[
        #align(horizon)[
     *Need to find a trade-off!* 
     ]
    ]
])
#slide(title: "Experiments: how to choose k?")[

    #image("images/k_choice_table.png")
    #pause
    *$k$ = 47*
    #pause
    - overrepresented fraction $lt$ 5%
    - rare fraction twice smaller
      
]
#slide(title: "Experiments: evaluation strategies")[
  
 Evaluated on two levels:
  
 + Computational efficiency:  
      - `AlfaPang` #emoji.swords `pggb` first two steps (`wfmash+seqwish`) 
 
 + Computational efficiency and output graph properties:
      - `AlfaPang+` = `AlfaPang + smoothxg + gfaffix` (`pggb` last two refinement steps) 
      
      - `AlfaPang+` #emoji.swords  `pggb` #emoji.swords `Minigraph-Cactus`
]

#slide(title: "Experiments: computational efficiency")[
  #image("images/performance_1.png")
  #image("images/performance_2.png")
]

#slide(title: "On the output graphs: graph topology")[
  *Goal:* measure the complexity of the produced pangenome graphs
  
  $=>$ comparison of the number of nodes and edges
  #columns()[
    #image("images/table_nodes.png")
    #image("images/table_edges.png")
  ]
  - We notice again that `MiniGraph-Cactus` performs worse: $17-45%\/10-58%$ _S. cerevisia_/_E. coli_ more *nodes* than `AlfaPang+`
  #image("images/table_nodes_length.png")
  - `AlfaPang+` shows an *higher rate compression*
]
#slide(title: "On the output graphs: graph similarity")[
      #image("images/table_aligned_pos.png")
      - number of aligned pairs: `AlfaPang+` larger than the others
      - A more precise analysis showed that:
        - `AlfaPang+` was able to find $96-99%$ of the pairs found by `pggb`
]
#focus-slide()[
  The end
  #pause
  ... almost
]
= Conclusion
#slide(title:"Summing up")[
  #linebreak()
  - We introduced `AlfaPang`: novel algorithm for building pangenome graphs
    - their structure is strictly defined by the $k$-completeness and $k$-faithfulness properties
    
  #pause
  
  - Runtime and memory usage scale *linearly* with the number of genomes
    #pause
    $=>$ process much larger sets than the other state-of-art alternatives
  #pause  
  - It is really *sensitive* to *sequence similarity*
    - deciding if given fragments should be aligned in a pangenome graph is somewhat arbitrary
 
]
#slide(title:"Summing up")[
   #pause  
  - `MiniGraph-Cactus`: differs in the assumptions on how the pangenome graph should look like
    - Make possible to reduce the number of sequence alignments 
    #pause
    $=>$ better computational efficiency than `pggb`
    #pause
    but still *worse* than `AlfaPang` in terms of *memory usage*
    
  - Close relationship with de Bruijn graphs
  #pause
  $=>$ Excessive *entanglement* in areas representing *low-complexity sequence regions*
  #pause
    - Thus removed by a refinement step by means of `smoothxg` tool #pause $=>$ `AlfaPang+`
    #pause  
        - This step dominates the computation time
  #pause $=>$ More precise tuning on parameters could allow to reduce it
]
#ending-slide(title: "The end")[Thanks for the attention!]
// #show: appendix

// = References

// #slide(title: "References")[
//   #bibliography("bibliography.bib", title: none)
// ]
