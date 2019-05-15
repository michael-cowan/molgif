![C60](gifs/C60.gif)

# molgif

create smooth gifs of rotating molecules

---

## Example

    import molgif
    import ase.build

    # load in molecule (ase.Atoms object)
    molecule = ase.build.molecule('biphenyl')

    # specify save path
    save_path = 'biphenyl.gif'

    # create rotating gif
    # auto_rotate: rotates molecule for better view
    molgif.rot_gif(molecule, save_path, auto_rotate=True)

![biphenyl](gifs/biphenyl.gif)

---

## Requirements

- ase
- matplotlib
