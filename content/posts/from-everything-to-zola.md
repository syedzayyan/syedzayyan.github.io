+++
title="From Everything to Zola"
date="2022-10-29T23:00:01.000Z"

[taxonomies] 
tags = ["thoughts"]
+++

I have a habit of rebuilding this blog often. Usually happens when I am trying to learn a new language or thing. This is like the comfort project that I can come back to. I used to think that this was a bad habit and just a thing I do when I am anxious. Surprisingly it's not. It's actually quite the opposite. David Thomas (the author of [Pragmatic Programmer](https://pragprog.com/titles/tpp20/the-pragmatic-programmer-20th-anniversary-edition/)) mentioned in his [blog](https://pragdave.me/blog/2020/05/05/advice-on-learning-a-language.html) that it's easier to learn a language that way.

The first iteration was hosted on Heroku, and built using Django-Rest. I built the front end using React and used MUI. It was something of a mess but it was easy to build in a few days. Everything went fine. I got to pick up React and its Hooks and work with a Rest API. It was during the first wave of COVID.

Quickly I learnt that React when it's not an SPA has problems with SEO. Since Google does not sweep JavaScript and only looks at the HTML it missed all the meta tags and title tags when using something like Helmet. I moved to NextJS. NextJS did not require a backend as such. I made all cool animations using framer motion. I added comments, users via NextAuth, Google OAuth and used AWS’s Amazon DocumentDB. I forget the name. I maintained this facade for a while but I didn't like the website as such as it still used MUI and didn't have my kind of style. It worked though so I kept it on and added sad cringe poetry. Yes I deleted the whole project from GitHub. But I learnt NextJS and [landed a freelancing gig](https://drreiter.com/). My only one till date.

I moved to Hugo because I was learning GoLang at the time and the themes were simple enough and just what I wanted. I know Hugo doesn't really need GoLang under the hood but emotions are illogical. I did write more though. Like good articles. I didn't know GitHub Actions, so I manually built the site deployed to GitHub Pages. This went on until I moved back to NextJS with a markdown converter and all on GitHub. I hosted the app on Vercel so all was taken care of.

And then enter now. Moving between countries, universities, job rejections later I have a lot of free time. I wanted to learn Typescript. Remember Pragmatic Programmer? Yeah it's my inspiration now for all things programming. I converted my [well plate app](wellplatelayout.xyz/) to react typescript from Svelte. Because I found this website seemed to daunting to rebuild again that too from NextJS.

After finishing that, I wanted to rebuild my site because I don't like Vercel as well. NextJS seems to locked to a Vendor. Enter Zola. Its simple. The lack of likeable themes meant I needed to learn the template which is actually pretty simple. Its on Zola now. I intend to publish a theme that I frankensteined from various themes for my site. Hopefully, everything goes alright contributing to open source for the first time.