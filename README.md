# Proton Disorder
## Overall approach
So I have a kind of complicated task ahead of me. The way I was approaching things before were very disorganized, and I was not keeping in mind how ice works. Before I was just randomly rotating the protons, which is not really what we want to do.

Instead, there are four main positions that the proton can be in. Or, there are 4 neighboring oxygen atoms, meaning each oxygen has 4 _links_ between them. Each link will contain just 1 proton. So I need an oxygen atom, that will contain a list of neighboring oxygens. I then need a list of all the links in the ice, with the oxygens that it contains. Then, put a hydrogen into each link and assign is to one of the oxygens. Do it randomly. Then, shake the links such that each oxygen in the ice has 2 protons assigned to each, and each link contains only 1 proton. Then calculate the dipole moment, and then shake again until you reach a dipole moment of 0.

That's the basic gist of what I need to do. As for an approach, I need to have an oxygen object that contains the number of neighboring oxygens, along with a number of protons associated with that oxygen. I then need a links object, which will contain the two oxygens that make up that link, and then the proton and which oxygen own the proton. I think. 

## Notes
I initially wanted to have the link object just be associated with the oxygen atoms. However, I think that I need a list of links, such that when I go to shake I can just iterate through the list rather than have to go through all the oxygen atoms and the associated link. 

I think that each link need to keep track of the oxygen atoms associated with it. I can do that through an index of the oxygen atom. So then I don't really need oxygens to know about the links, but I do need them to know how many links they have associated to them. I can just use a counter in this case, something like `nneighbours`. I think that should work nicely.

### init_hydrogens
So I've had a bit of a speed bump. In this function, we are trying to place a hydrogen in each link, and have that hydrogen be associated to one of the oxygens within that link. The way that this is done in the original code is with these objects called half links, where each oxygen has a half link and they contain information like how many bonds that oxygen had and whether the proton within that bond is associated with it. It's just confusing why this was done. I'm trying to come up with a way to sort of simplify this set up, or at least reduce the number of objects being use.

_Update 6/15/23_

I played around with it a little more, and I decided to not use a `HalfLinks` class. Instead, my Links class tracks both oxygens, and then whether the bond is on one oxygen or the other. In addition, each oxygen tracks each of the links that it is associated with via index. So there are two main lists that I have - `ice` and `links`. Ice contains all of the oxygen objects, links contains all the `link` objects. All the `links` objects contain indexes that relate to the `ice` array, and vice versa. This is my workaround for not having pointers in python. This whole set up would be a lot easier to implement if I had pointers for everything.

### shake_bonds
This method was relatively easy to implement. It helped that once I had set up this relation between the two arrays, I could update the two quite readily. So yeah - it does what it needs to do. 

### get_dipole
I'm a little stumped on this, but that is okay. My oxygens need to know the location of the hydrogens associated with it. However, the oxygen objects that I have do not know anything about the hydrogens. All they know are neighbours, through link objects. So, what I can do is parse through the links, and see which ones have this oxygen as having a bond. I then need to get the coordinates of the neighbouring oxygen, because that will allow me to get the direction of the hydrogen. Once I have the direction for both 

# Attempt v2
So I took a moment, stepped away and decided to approach things a little differently. I was making individual oxygen and links objects, and having them point to each other via indexes in arrays. But I think the fact that I'm using objects is making things more complicated than it needs to be. So instead, I'm going to use a 2D array (which granted doesn't seem like it would be less complicated). The reasoning here is that it should allow me to to access things easier, without the added complexity of classes.

## Array structure
Here is the main layout of the 2D array:
```
[
    [[x,y,z],[x,y,z],...], # Coordintes for oxygen atoms
    [
        {1:True,2:False,...},
        {0:False,2:True,...},
    ] # The links for each oxygen
] 
```

The first part of the array is an array of 3-tuples that contain the oxygen coordinates. This is parsed in by the xyz/csv file.

The second part is really where the magic happens. The second array is filled with dictionaries, each of which has four entries for each link associated with oxygen. The original idea was to use tuples, but since tuples are immutable in python, that makes it quite difficult to operate on.

One problem - instantiating hydrogens. The link objects are related to each other. So I need to make sure that I can index them randomly, but that if I change one I must change the other.