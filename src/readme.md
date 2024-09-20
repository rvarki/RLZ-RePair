Tools from the repair library designed to avoid the generation of very deep rules

Both repair versions take as an additional parameter an integer providing the amount of available RAM in MB. The use of such limit only affects the running time, not the shape of the grammar


* **irepair**, **idespair** RePair compression and decompression for sequences of 32 bit integers. **irepair** is used by BigRePair; **idespair** should be equivalent to **repair/idespair**.


* **repair**, **despair** RePair compression and decompression for sequences of characters using the original RePair format (including an alphabet map, see despair.c)

