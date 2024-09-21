+++
title="Rust is for lobotomised infants"
date="2024-09-23T23:00:01.000Z"

[taxonomies] 
tags = ["thoughts", "rants", "programming"]
+++

### Introduction

Boom Clickbaited! Hold up your pitchforks, Rust folks. C folks I will disappoint you too, just like I did my parents—so hold off on the party poppers. Rust is quite a weird language coming from the world from Python and even weirder from the world of C. I personally come from the world Javascript and probably some Python (more Pytorch/ML than idiomatic Python but you get the point). The type system in Javascript is [Wat?](https://www.destroyallsoftware.com/talks/wat) and for Python it's fairly non-existent (I know but like ... you get the point). Rust is too far right into types. Recently I worked with implementing a server using the Postgres Wire Protocol and used Rust for it. This blog is my rant, realisation and probably some note to self so that I can laugh at myself a few years down the line.

### Memory and Pointers

Firstly, Rust considers raw pointers to be "unsafe". C folks are probably rolling in laughter, but not everyone likes to do mental gymnastics with memory. Before you script kiddies shit yourselves at the mention of memory, here's a short summary (not that I am not a script kiddy myself). Everything in the computer is memory. Like reading this boring article you are on is also on memory (in some form) that the browser is reading from. Programming is basically moving and handling this memory. C and Assembly gives you access to this memory and it is up to you allocate some for yourself, and free it later on. In assembly you work with memory addresses. You can think of memory as like dairy milk bars. Each block is a memory cell, and you store your stuff there. The 'address' is like a label that tells you where to find your stuff later.

In C, you get a lovely abstraction called pointers. Pointers are variables that hold the memory addresses. So you can just pass around this pointer instead of remembering addresses and do crazy mental gymnastics. Yes, I know in those days people did exactly that, but we have evolutioned into the age of abstractions!

### Thinking of memory in a way to avoid gymnastics

A good way to think about memory is in the form of ownership. I don't know how common it is but this [blog](https://www.chiark.greenend.org.uk/~sgtatham/cdescent/) mentions it so I am hoping it is common. Basically if we have same pointer used by various parts of the program, we consider one part/function the owner of this pointer. The owner has the responsibility of free-ing this pointer and the memory in that pointer. You don't want to free the memory space too early as it may lead to segfaults as other parts accesssing this memory may never find it. Rust takes this concept and builds a language around it. Which is pretty cool for people like me who can hurt themselves sitting down on a chair. Rust has this ownership model which is better demonstrated in code:

```Rust
fn main() {
    let r; // Declare a reference (but it's uninitialized)
    {
        let x = 5;  // x is created here
        r = &x;     // Borrow x and assign the reference to r
        // r can be used here within the inner scope
    }
    // x goes out of scope here
    // r is now a dangling reference because x is no longer valid
    println!("r: {}", r); // ERROR: `x` does not live long enough
}
```

```Rust
error[E0597]: `x` does not live long enough
 --> src/main.rs:6:13
  |
4 |         let x = 5;
  |             - `x` dropped here while still borrowed
5 |         r = &x;
6 |     }
  |     - borrow later used here

```

```Rust
fn main() {
    let x = 5;   // x is declared in the outer scope
    let r = &x;  // r borrows x
    println!("r: {}", r); // Now r is valid because x is still in scope
}
```

When you force someone to follow a paradigm, you have a problem! Now raw pointers are not really dangerous, as people would like to say (unless you are a newb like me). The computer works with memory, and a programmer's sole responsibility is to write good software by using memory in an efficient manner. There's obviously hundred different ways, you could be clever with your memory. From garbage collection methods where you just let your memories dangle and run algorithms periodically to clean up unused memory. And all of them are clever. Rust re-introduces this ownership model.

### Rust claims I don't necessarily agree with

Rust touts itself as a low-ish level language, but wants to behave like Java while the language expression feels like a functional one. It's a sort of hybrid new identity that feels utterly wrong from the decades of language that it precedes. I am not saying its bad by any means, I am just saying its different. Although I do take some issue to some of the things Rustaceans claim. Like Rust being blazingly fast. With that much safety checks, it just can't be blazingly fast as C. I know LLVM optimisations exist but let me have my opinion in peace. 

Additionally, achieving fast in Rust is doing mental gymnastics with ownership and borrowing. And then we fall into the same pitfalls as C. So its not really fast and its not really safe, because most pesky time consuming bugs end up being logical errors and I can't switch to another brain. I have not mastered Rust or C or programming for that matter. Maybe somewhere down the line I'll have an opinion on the sorts of paradigm I like.

### Then what is Rust really good at then? 

Make you think differently. Its never really a waste of time to learn new languages (probably not Javascript) and if its one that shifts your thinking completely its always worth a shot. I guess I just want to say programming is mental gymnastics. While language like Python make your life easy by abstracting away all of memory management, it's still your responsibility to write good software. So, the question becomes: Do you want to do bare-metal memory gymnastics with C? Or would you rather take on memory gymnastics within Rust’s set paradigm—where you’re guided (sometimes forced) into safety but still have the power to go "unsafe" when needed? Maybe it just depends. 

So here I am, learning C, because sometimes you need to embrace the jank and learn to put on your big boy pants up.