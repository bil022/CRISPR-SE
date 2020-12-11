
# Genome editing tool upgrade:

- Replacing CRISPOR search engine with CRSIPR-SE:

CRISPOR use BWA for mismatch search with following parameters:
```
bwa aln -o 0 -m 1980000 -n 4 -k 4 -N -l 20 EcoliE84.fa $input.fa > $input.sa
bwa samse -n 60000 EcoliE84.fa $input.sa $input.fa 
```
## The following scripts use CRISPR-SE with same output format as BWA: 
```
se --index -sr $input
se --build -m 5 -v -r EcoliE84 -q $input | crispor-se.pl --format | sort | crispor-se.pl --parse --ref EcoliE84 --query $input 
```
## For example:
For the query of GATGGCGTTTAATCGCCTTCCGG with Ecoli from Chromosome:254218-254240, reverse strand:

CRISPOR with [CRISPR-SE](http://renlab.sdsc.edu/CRISPR-SE/crispor/crispor.py) reported two off-targets, where [CRISPOR](http://crispor.tefor.net/crispor.py) can only report the first one:
```
guide:      GATGGCGTTTAATCGCCTTC CGG
off-target: GATGGCCATGAATGGCCTTC AGG
                  ** *   *      
CFD Off-target score: 0.000000
MIT Off-target score: 0.13
Position: Chromosome:782981-783003:+
Distance from target: 0.529 Mbp
```
```
guide:      GATGGCGTTTAATCGCCTTC CGG
off-target: GAAGGCGATTAAACGCCATC CGG
              *    *    *    *  
CFD Off-target score: 0.263736
MIT Off-target score: 0.12
Position: Chromosome:254221-254243:+
Distance from target: 0.000 Mbp
```


Thanks for visiting [The Markdown Guide](https://www.markdownguide.org)!

This Markdown cheat sheet provides a quick overview of all the Markdown syntax elements. It can’t cover every edge case, so if you need more information about any of these elements, refer to the reference guides for [basic syntax](https://www.markdownguide.org/basic-syntax) and [extended syntax](https://www.markdownguide.org/extended-syntax).

## Basic Syntax

These are the elements outlined in John Gruber’s original design document. All Markdown applications support these elements.

### Heading

# H1
## H2
### H3

### Bold

**bold text**

### Italic

*italicized text*

### Blockquote

> blockquote

### Ordered List

1. First item
2. Second item
3. Third item

### Unordered List

- First item
- Second item
- Third item

### Code

`code`

### Horizontal Rule

---

### Link

[title](https://www.example.com)

### Image

![alt text](image.jpg)

## Extended Syntax

These elements extend the basic syntax by adding additional features. Not all Markdown applications support these elements.

### Table

| Syntax | Description |
| ----------- | ----------- |
| Header | Title |
| Paragraph | Text |

### Fenced Code Block

```
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```

### Footnote

Here's a sentence with a footnote. [^1]

[^1]: This is the footnote.

### Heading ID

### My Great Heading {#custom-id}

### Definition List

term
: definition

### Strikethrough

~~The world is flat.~~

### Task List

- [x] Write the press release
- [ ] Update the website
- [ ] Contact the media
