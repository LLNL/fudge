.. _glossary:

********
Glossary
********

.. if you add new entries, keep the alphabetical sorting!

.. glossary::

    component
        A GNDS 'reaction' object contains several types of physical data, including the crossSection, the outputChannel
        with a corresponding Q-value and list of products, and for each product a multiplicity and distribution.
        Each of these data types is a 'component' which contains one or more data 'form' (see below). Each form
        in the component should have a label corresponding to a 'style' defined in the reactionSuite styles section.

    evaluation
        A complete set of parameters (e.g. the cross sections) defined at all energies, even when no experimentally measured data are available.
        It is the task of the evaluator to assess the most probable value of a parameter at any energy, resolving issues of discrepant measurement,
        assigning values (by an educated guess or based on model calculations) where no data are available and providing data in computer-readable format.

    form
        A form is one representation of a type of data. Forms are collected inside of a component. For example, the crossSection component
        may contain an evaluated form (which might be 'resonancesWithBackground') along with several processed forms including
        'XYs1d' (after resonance reconstruction and/or Doppler broadening) and 'gridded1d' (after grouping).
        Each form should have a label corresponding to a 'style' defined in the reactionSuite styles section.

    PoPs
        Properties of Particles (PoPs) is a database for storing particle information such as mass, spin, halflife, etc.
        A PoPs database may be stored inside of a reactionSuite or may be stored as a separate stand-alone file.

    protare
        Contraction for 'PROjectile + TARget + Evaluation'. Also known as a reactionSuite, a protare contains all possible reactions involving this projectile/target combination.

    reactionSuite
        Contains all possible reactions involving a projectile/target combination. May contain additional information such as summed reactions. Also sometimes called a Protare

    style
        Defines a type of data stored inside a GNDS file. The most common style is 'evaluated' (corresponding to data added by an evaluator).
        Other styles include 'crossSectionReconstructed', 'heated', 'griddedCrossSection', etc. Each style defines some details about the data + how
        they were derived, and each style has a unique label used to associate that style with data elsewhere in the file.
        Also see the definitions of 'component' and 'form'.

    suite
        A suite in GNDS contains a collection of similar objects. For example, the 'reactions' suite contains a list of 'reaction' instances
        while a 'nuclides' suite (in the PoPs particle database) contains a list of 'nuclide' objects.
