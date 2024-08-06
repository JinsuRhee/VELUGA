"""
This file is automatically generated from IDL files.
"""

def r_gal():
    """
    ; Load Galaxy/Halo Catalog Data.
    ; This method retrieves galaxy or halo catalog data for a given snapshot and object ID.
    ;
    ; Parameters
    ; ----------
    ; snap0 : int
    ;	Snapshot number
    ;
    ; id0 : int
    ;	Object ID. Use a negative value to retrieve all objects in the snapshot.
    ;
    ; horg: 'h' or 'g'
    ;	A flag to specify the object type. Galaxy for 'g' and Halo for 'h'
    ;	Default is 'g'
    ;
    ;
    ; Returns
    ; -------
    ; Structured_array
    ;	A structured array containing information about the objects.
    ; g = veluga->r_gal(100, 1)
    ;
    ; Examples
    ; --------
    ; IDL> g = veluga->r_gal(100, 1)	; Read the galaxy with ID=1 at the snapshot of 100
    ; IDL> PRINT, g[0].ID 				; Print its ID
    ;
    ; IDL> h = veluga->r_gal(200, -1, horg='h') ; Read all halos at the snapshot of 200
    ; IDL> PRINT, h[0].mvir 					; Print the virial mass of the first halo
    """
    pass

