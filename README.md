# Proton Disorder
## Overall approach
So I have a kind of complicated task ahead of me. The way I was approaching things before were very disorganized, and I was not keeping in mind how ice works. Before I was just randomly rotating the protons, which is not really what we want to do.

Instead, there are four main positions that the proton can be in. Or, there are 4 neighboring oxygen atoms, meaning each oxygen has 4 _links_ between them. Each link will contain just 1 proton. So I need an oxygen atom, that will contain a list of neighboring oxygens. I then need a list of all the links in the ice, with the oxygens that it contains. Then, put a hydrogen into each link and assign is to one of the oxygens. Do it randomly. Then, shake the links such that each oxygen in the ice has 2 protons assigned to each, and each link contains only 1 proton. Then calculate the dipole moment, and then shake again until you reach a dipole moment of 0.

That's the basic gist of what I need to do. As for an approach, I need to have an oxygen object that contains the number of neighboring oxygens, along with a number of protons associated with that oxygen. I then need a links object, which will contain the two oxygens that make up that link, and then the proton and which oxygen own the proton. I think. 

## Notes
I initially wanted to have the link object just be associated with the oxygen atoms. However, I think that I need a list of links, such that when I go to shake I can just iterate through the list rather than have to go through all the oxygen atoms and the associated link. 

I think that each link need to keep track of the oxygen atoms associated with it. I can do that through an index of the oxygen atom. So then I don't really need oxygens to know about the links, but I do need them to know how many links they have associated to them. I can just use a counter in this case, something like `nneighbours`. I think that should work nicely.