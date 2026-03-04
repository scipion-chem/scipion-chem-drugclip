=======================
Biofold plugin
=======================

**Documentation under development, sorry for the inconvenience**

This is a **Scipion** plugin that offers tools from the
`Boltz <https://github.com/jwohlwend/boltz>`_, and `Chai <https://github.com/chaidiscovery/chai-lab>`_, as well as imports for `AlphaFold3 <https://github.com/google-deepmind/alphafold3>`_, `Protenix <https://github.com/bytedance/Protenix>`_, Boltz and Chai server results.
These tools will make it possible to carry out different functions to fold proteins.


==========================
Install this plugin
==========================

You will need to first install
`Scipion3 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_  and
`Scipion-chem <https://github.com/scipion-chem/scipion-chem>`_ to run these protocols.


1. **Install the plugin in Scipion**

Biofold is installed automatically by scipion.

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or

.. code-block::

    scipion3 installp -p scipion-chem-biofold


- **Developer's version**

    1. **Download repository**:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-biofold.git

    2. **Switch to the desired branch** (master or devel):

    Scipion-chem-biofold is constantly under development and including new features.
    If you want a relatively older an more stable version, use master branch (default).
    If you want the latest changes and developments, user devel branch.

    .. code-block::

                cd scipion-chem-biofold
                git checkout devel

    3. **Install**:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-biofold --devel




