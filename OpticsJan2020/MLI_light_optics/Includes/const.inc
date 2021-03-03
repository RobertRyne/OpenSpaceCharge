cryne 7/7/2002
cryne this common block is used to store strings that have constant
cryne values assigned to them.
cryne nconst is equal to the number that have been assigned
cryne This info should really be included in elments.inc, but I am
cryne keeping it separate to avoid changes to the existing MaryLie
      parameter (nconmax=1000)
      character*16 constr
      common/conarray/constr(nconmax),conval(nconmax),nconst
