# EC Index Calculus Benchmarks
This is the code used to create instances derived from a Weil Descent on the fourth Summation polynomial for an Index Calculus attack on binary elliptic curves.

There are already some benchmarks in the folder named 'benchmarks' which were created using
```bash
sage find_points.sage 15 5 10 > launch.sh
chmod +x launch.sh
./launch.sh

sage find_points.sage 17 6 10 > launch.sh
chmod +x launch.sh
./launch.sh

sage find_points.sage 19 6 10 > launch.sh
chmod +x launch.sh
./launch.sh
```

You can create other benchmarks in a similar way, by choosing the arguments for the ```find_points.sage``` script, which correspond to parameters ```n```, ```l``` and ```nb```. This will create a total of ```2⋅nb``` instances where at least ```nb``` instances will have at least one solution. 


