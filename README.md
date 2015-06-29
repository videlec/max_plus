# A small Python library to deal with max-plus matrices relations

To run the library you need to install Sage (http://sagemath.org). Then
clone this repository (or download all the files). Go to the directory
where you downloaded the file. Start Sage and run

    sage: %runfile int_max_plus.pyx
    sage: %runfile max_plus.py

Then you can for example list of relations for the 2x2 triangular matrices

    sage: relations_tri(2)
    xyyxxyxyyx = xyyxyxxyyx
    xyyxxyyxxy = xyyxyxyxxy
    xyyxyxxyyx = xyyxxyxyyx
    xyyxyxyxxy = xyyxxyyxxy
    xxxxyyyyxxyyxy = xxxxyyyyxyxyxy
    xxxxyyyyxyxyxy = xxxxyyyyxxyyxy
