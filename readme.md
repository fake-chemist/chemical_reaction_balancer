# Chemical Equation Balancer 

This tool is meant for balancing simple chemical equations

## Motivation 

Chemistry students everywhere in Highschool and College alike are learning to 
balance equations. The online tools I've found have adds and I don't like that 
so I made a quick tool for anyone to use.  


## Usage

```bash 
python reaction_balancer.py "C2H5OH + O2 -> CO2 + H2O"
```

## Output

```bash
C2H5OH + O2 -> CO2 + H2O
_______________________________
C2H5OH + 3O2 -> 2CO2 + 3H2O
```

### Dependencies 
  -pulp   
  -pytest
