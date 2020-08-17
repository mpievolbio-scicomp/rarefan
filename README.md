# RepinPop

## Create the conda environment

```
$> conda env create --file=environment.yml
```


## Build the java code:
`RepinPop` requires at least java version 11. Building is done by
[`gradle`](https://gradle.org).

```
$> cd REPIN_ecology/REPIN_ecology
$> gradle build
```

## Launch the server
Our application is served as a web form where users can upload their sequence
files, set the parameters, run the analysis and visualize the results. The web
form is implemented in a jupyter notebook and served with the `voila` framework.
To launch the server, run

```
$> voila --enable_nbextensions --port=<YOUR_FREE_PORT> --no-browser --Voila.ip=0.0.0.0
```

NOTE: You have to specify the port on which to run the service. The last
argument will make the server listen to requests coming from any IP. Restrict
the range of IPs or leave out this argument if this behaviour is unwanted. For
more advanced server configurations see the [`voila`
documentation](https://voila.readthedocs.io/en/stable/deploy.html).




