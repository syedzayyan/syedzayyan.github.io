+++
title="Driving into a tree, in a Rusty car"
date="2024-07-08T08:15:51.000Z"

[taxonomies]
tags=["rants", "rust", "code", "discussion"]
+++

I like Rust, at least I want to on most days, but implementing trees in Rust is interesting. First of all, why do I even need a tree. I am writing a cheminformatics toolkit and am implementing a SMARTS parser. I am going the ubiquitous route of building the abstract syntax tree (AST) to match a molecule structure. Easy! (Basically I am stealing [OpenBabel's SMARTs Parser ](https://github.com/openbabel/openbabel/blob/master/include/openbabel/parsmart.h)because I am an impostor amongst real programmer) But, Rust's strict ownership rules make it really difficult to have tree structures implemented in the way it would be written in a better language like C/C++. I will shuffle between hating and loving Rust so hold your pitchforks till the end. So, in C++ a `BondExpr` will be defined as something like:

```C++
  typedef union _BondExpr {
    int type;
    struct
    {
      int type;
      union _BondExpr *arg;
    }
      mon;
    struct
    {
      int type;
      union _BondExpr *lft;
      union _BondExpr *rgt;
    }
      bin;
  } BondExpr;
```

Unions don't exist in Rust, I mean they do but you need `unsafe`, and that point might as well use C and shoot yourself in the foot. So we can use the better thing: `enum`! 

```Rust
pub enum BondExpr {
	Monadic {
		bond_type: BondTypes,	
		arg: BondExpr,
	},
	Binary {
		bond_type: BondTypes,
		left: BondExpr,
		right: BondExpr,
	},
}
```

Rust being Rust, you can't have that. We get an error:

`recursive type BondExpr has infinite size`

I mean that makes sense but how do I just point to a variable that's allocated somewhere and I promise I'll clean it. You don't. In Rust you get a `Box`, because it thinks I can't take care of myself. I mean sure, I can't, but show some faith. Using boxes, we reach something that looks ugly, but something that works:

```Rust
pub enum BondExpr {
    Monadic {
        bond_type: BondTypes,
        arg: Option<Box<BondExpr>>,
    },
    Binary {
        bond_type: BondTypes,
        left: Option<Box<BondExpr>>,
        right: Option<Box<BondExpr>>,
    },
}
```

Now you might go, but why can't you just use a flat vector and store indices. Sure, great idea, except I am doing that already. 

```Rust
pub enum NodeData {
    Atom(AtomSpec),
    Bond(BondSpec),
    Unknown,
}
pub struct TreeNode {
    pub op_code: OpCode,
    pub data: NodeData,
    pub parent: usize,
    pub visit: bool,
}
pub struct SmartsPattern {
    pub nodes: Vec<TreeNode>,
    pub root: usize,
    pub smarts_string: String,
    pub chirality: bool
}
```

But I don't want to make a bond node if I can get away with just pointing to a `BondExpr`. And then you reach the problem where Bonds have a source and destination and you can often infer that the next node after the bond will be an atom node, but if you have an expression hogging up your clean vector, now you have the problem of going back fix your destination in the `BondSpec`. Now we could have a second vector for storing these second layer of expressions, but then when setting the indices it becomes a nightmare to read the code. Sure we can make our own types and do a casting but eh.

So yeah. Boxes and Options it is with a bunch of unwraps when matching the molecule. I am not complaining, it's just what it is in Rust, or if you have a better idea, please by an mean contribute to my [molrus](https://github.com/syedzayyan/molrus).
