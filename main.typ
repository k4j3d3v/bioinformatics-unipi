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
  progress-bar: false,
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
  #image("images/download.png")
]
#slide(repeat: 3, self => [
  #let (uncover, only, alternatives) = utils.methods(self)

  At subslide #self.subslide, we can

  use #uncover("2-")[`#uncover` function] for reserving space,

  use #only("2-")[`#only` function] for not reserving space,

  #alternatives[call `#only` multiple times \u{2717}][use `#alternatives` function #sym.checkmark] for choosing one of the alternatives.
])

#slide(title: "Theorems")[
  Theorems can be created with the `#theorem` command. Similarly, there are `#proof`, `#definition`, `#example`, `#lemma`, and `#corollary`. \
  For example, here is a theorem:
  #theorem(title: "Important one")[
    Using theorems is easy.
  ]
  #proof[
    This was very easy, wasn't it?
  ]
  A definition already given by well-known mathematicians @Author1978definition is:
  #definition(title: "Important stuff")[
    _Important stuff_ is defined as the stuff that is important to me:
    $
      exp(upright(i) pi) + 1 = 0.
    $
  ]
]

#slide(title: "Equations")[
  Equations with a label with a label will be numbered automatically:
  $
    integral_0^oo exp(-x^2) dif x = pi/2
  $<eq:important>
  We can then refer to this equation as @eq:important.
  Equations without a label will not be numbered:
  $
    sum_(n=1)^oo 1/n^2 = pi^2/6
  $
  Inline math equations will not break across lines, which can be seen here: $a x^2 + b x + c = 0 => x_(1,2) = (-b plus.minus sqrt(b^2 - 4 a c))/(2 a)$
]

#show: appendix

= References

#slide(title: "References")[
  #bibliography("bibliography.bib", title: none)
]
