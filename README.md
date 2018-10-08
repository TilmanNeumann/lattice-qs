# lattice-qs

A study of a "lattice quadratic sieve" proposed by R.D.Silverman in https://www.mersenneforum.org/showthread.php?t=14080.

The algorithm is not competitive with SIQS, rather comparable to a good CFrac implementation. The biggest number I factored yet
using it had 220 bit, and that took around 40 minutes.

Nonetheless, eventually the implementation of the lattice sieve might be interesting to some people.
But beware: I wrote it without having seen Pollard's "The lattice sieve" paper, so it may look quite different (and probably worse) than usual approaches. My best source of information was [Franke, Kleinjung 2005: CONTINUED FRACTIONS AND LATTICE SIEVING].

## Releases

* v0.1 the one and only?


## Getting Started

Clone the repository, create a plain Java project importing it, make sure that 'src' is the source folder of your project, and add the jars from the lib-folder to your classpath. 

You will need Java 8 or higher for the project to compile.

There is no documentation and no support, so you should be ready to start exploring the source code.


## Authors

Tilman Neumann


## License

This project is licensed under the GPL 3 License - see the [LICENSE](LICENSE) file for details


## Credits

www.mersenneforum.org

