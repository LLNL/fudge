
Introduction
============

Properties of Particles (or PoPs) is a database for storing particles and their reaction-independent properties.
The database is meant to be general enough to describe all relevant particles and their properties, including mass,
charge, spin, parity, half-life and decay properties. It stores not only fundamental particles,
nucleons and excited nuclear states, but also excited *atomic* states that can de-excite through atomic relaxation.

Every particle stored in a PoPs database requires a unique string id, that can be used to refer to or look up the
particle. Particle ids should follow naming conventions to make them easier to
use and understand.

PoPs works together with the Generalized Nuclear Database Structure (GNDS). In general, PoPs stores information about
individual particles, while GNDS stores information about interactions between those particles.
