let scaffold = systems[0].strands.filter(s=>s.getLength()>500)[0];
let domains = []; 
let domain = [];

scaffold.getMonomers().forEach(base => {    
    let pair = base.pair;
    // if we have a pair and domain is not empty 
    if(pair){
        // we check that pair has the same strand as the bases in domain
        if(domain.length == 0 || domain[0].strand === pair.strand){ 
            domain.push(pair);
        } else {
            //start a new domain if we have a pair with a different strand
            if(domain.length > 0){
                domains.push(domain);
            }
            domain = [pair];
        }

    }
});

// now we want to go over the domains type is where base info is stored as a string  
// so we can convert domains to sequences 
// because of 5' to 3' direction we need to reverse the sequences
let sequences = domains.map(domain => domain.map(base => base.type).reverse().join(''));

//now lelt's print the sequences 1 at a line 
makeTextFile("domains.txt",  sequences.join('\n'));