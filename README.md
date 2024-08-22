# Decoding CSIDH: A Comprehensive Guide to Isogeny-Based Cryptography
This is the uncompiled version of my Master Thesis at Leiden University.
(The compiled version can be viewed here as the file `build_latex/main.pdf`)

Also, a huge thanks to [VimTeX](https://github.com/lervag/vimtex) for making it so easy to compile my tex files and interacting with [zathura](https://github.com/pwmt/zathura) enabling hot reload. It allowed me to edit my thesis in [NeoVim](https://github.com/neovim/neovim) whilst only needing to save the file in order to see the changes in my pdf-viewer.

Some notes:
1.  I am a huge fan of vectorised pictures, this allows my thesis of 67 pages to be 587 KiB.
2.  I designed the file structure to allow me to save copies of the current version without copying all plots and pdfs. I did this by continuously working in the `current` folder, where I can generate all the plots I use as well as compile the thesis. The output files would be saved outside of the `current` folder. This allowed me to save copies and backups of the `current` state of my thesis by copying the `current` folder, only containing the files that generate the final version.
3.  Before compiling my thesis, one can generate all plots by running `cd current` and `python generate_plots.py`.
4.  To compile my thesis I used the compiler `xetex` (it completely compiles without any errors and warnings) after installing [`texlive-most`](https://archlinux.org/groups/x86_64/texlive-most/) on Arch Linux. Specifically, I used the VimTeX settings:
    ```
    let g:vimtex_view_method = 'zathura'

    let g:vimtex_compiler_method = 'latexmk'
    let g:vimtex_compiler_latexmk_engines = {
    \   '_' : '-xelatex',
    \}
    let g:vimtex_compiler_latexmk = {
    \   'aux_dir' : '../build_latex',
    \   'out_dir' : '../build_latex',
    \}
    ```
    to compile the file in `current/main.tex` to `build_latex/main.tex`.

## Sage Code
In the file `current/display_code.sage` I have added the code to compute the example in Section 5.2 of my thesis.

Upon removing every occurrence of `.montgomery_model()` in the code, one will get the example in Section 5.3.0 of my thesis.

The code was written with the mindset of providing an idea of how the theory translates into practice and isn't optimised for performance.
One can freely tweak the parameters in the code, like the prime `p` that defines the finite field you are working over as well as the ideals that are applied.
