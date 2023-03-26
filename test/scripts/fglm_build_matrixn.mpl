with(Groebner):

is_staircase_grevlex_generic:= proc (ns,rv,lmGb,vars)
description "Check if the staircase ns is DRL-generic for the last variable "
    "of vars. Returns true if it is and false if it is not.":
local xn,ns_sorted,lmGb_sorted, mon_ord, ind_ns,ind_lmGb,m,i,j,cmp:
    xn     := vars[-1]:
    mon_ord:= tdeg (op(vars)):
    cmp    := (a,b)->Groebner:-TestOrder (a,b,mon_ord):
    ns_sorted:= sort (ns,(a,b)->cmp (a,b)):
    lmGb_sorted:= sort (lmGb,(a,b)->cmp (a,b)):
    j:= 1:
    for i from 1 to nops (ns_sorted) do
        m:= ns[i]*xn:
        if (type (rv[m],indexed)) then
            if (m <> lmGb[j]) then
                printf ("%a goes out the staircase and does not land on "
                        "a leading monomial\n",m):
                return false;
            else j:=j+1:
            end if:
        end if:
    end do:
    return true:
end proc:

ToMSolve:=proc(F, char, vars, fname)
local i, fd, F2, str;
    fd:=fopen(fname, WRITE):
    for i from 1 to nops(vars)-1 do
        fprintf(fd, "%a, ", vars[i]):
    end;
    fprintf(fd, "%a ", vars[nops(vars)]):
    fprintf(fd,"\n");
    fprintf(fd,"%d\n", char);
    if char = 0 then
        F2:=map(f->sort(expand(numer(f)), order=tdeg(op(vars))), F):
        F2:=remove(member, F2, [0]):
        for i from 1 to nops(F2)-1 do
            fprintf(fd, "%a,\n", F2[i]):
        od:
        fprintf(fd, "%a\n", F2[nops(F2)]):
    else
        F2:=map(f->sort(expand(numer(f)), order=tdeg(op(vars))) mod char, F):
        F2:=remove(member, F2, [0]):
        for i from 1 to nops(F2)-1 do
            fprintf(fd, "%a,\n", F2[i]):
        od:
        fprintf(fd, "%a\n", F2[nops(F2)]):
    fi:
    fclose(fd):
#   str := cat("sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\\\n//g' ", fname):
    str := cat("sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\\\\\\n//g' ", fname):
    system(str):
end proc:

main:= proc (F,vars,char,str)
local G,lmGb,lmGbxn,ord,ns,rv,fd,b,i,i2,j,nr_col,nr_row:
local dense_idx,triv_idx,triv_pos,Mxn,dst,cpt;
    ord  := tdeg(op(vars)):
    G    := Groebner:-Basis (F,ord,characteristic=char):
    lmGb := map (g->LeadingMonomial (g,ord), G):
    lmGbxn:= select (g->degree (g,vars[-1])>0,lmGb):
    nr_row:= nops(lmGbxn):
    ns,rv:= NormalSet (G,ord):
    nr_col:= nops(ns):
    triv_idx:= [seq (i,i=1..nr_col)]:
    dense_idx,triv_idx:= selectremove (l->type(rv[ns[l]*vars[-1]],indexed),
                                       triv_idx):
    triv_pos:= map (l->rv[ns[l]*vars[-1]],triv_idx):
    dst:= map (l->0,dense_idx):
    ToMSolve (F,char,vars,cat("../input_files/",str,".ms")):
    fd:= fopen (cat("../fglm/build_matrixn_",str,".c"),WRITE):
    cpt:
    fprintf (fd,
             "#include \"../../src/msolve/libmsolve.c\"\n\n"
             "int main(void)\n{\n"
             "    int32_t la_option             = 2; // by default\n"
             "    int32_t use_signatures        = 0;\n"
             "    int32_t nr_threads            = 1;\n"
             "    int32_t info_level            = 0;\n"
             "    int32_t initial_hts           = 17;\n"
             "    int32_t max_pairs             = 0;\n"
             "    int32_t elim_block_len        = 0;\n"
             "    int32_t update_ht             = 0;\n"
             "    int32_t generate_pbm          = 0;\n"
             "    int32_t reduce_gb             = 1;\n"
             "    int32_t print_gb              = 0;\n"
             "    int32_t genericity_handling   = 2;\n"
             "    int32_t saturate              = 0;\n"
             "    int32_t colon                 = 0;\n"
             "    int32_t normal_form           = 0;\n"
             "    int32_t normal_form_matrix    = 0;\n"
             "    int32_t is_gb                 = 0;\n"
             "    int32_t get_param             = 0;\n"
             "    int32_t precision             = 128;\n"
             "    int32_t refine                = 0; /* not used at the moment */\n"
             "    int32_t isolate               = 0; /* not used at the moment */\n"
             "    files_gb *files = malloc(sizeof(files_gb));\n"
             "    files->in_file  = \"test/input_files/%s.ms\";\n"
             "    files->out_file  = NULL;\n"
             "    FILE *fh  = fopen(files->in_file, \"r\");\n\n"

             "    //  int32_t mon_order   = 0;\n"
             "    int32_t nr_vars     = 0;\n"
             "    int32_t field_char  = 9001;\n"
             "    int32_t nr_gens     = 0;\n"
             "    data_gens_ff_t *gens = allocate_data_gens();\n\n"
             "    get_data_from_file(files->in_file, &nr_vars, &field_char,&nr_gens,gens);\n\n"
             "    if (nr_vars != %a) return 101;\n"
             "    if (field_char != %a) return 102;\n"
             "    if (nr_gens != %a) return 103;\n\n"
             "    gens->rand_linear           = 0;\n"
             "    gens->random_linear_form = malloc(sizeof(int32_t)*(nr_vars));\n\n"
             "    param_t *param  = NULL;\n"
             "    mpz_param_t mpz_param;\n"
             "    mpz_param_init(mpz_param);\n\n"
             "    long nb_real_roots      = 0;\n"
             "    interval *real_roots    = NULL;\n"
             "    real_point_t *real_pts  = NULL;\n\n"
             "    param_t ** paramp = &param;\n"
             "    mpz_param_t *mpz_paramp = &mpz_param;\n"
             "    long *nb_real_roots_ptr = &nb_real_roots;\n"
             "    interval **real_roots_ptr = &real_roots;\n"
             "    real_point_t **real_pts_ptr = &real_pts;\n\n"
             "    int32_t *bld    = NULL;\n"
             "    int32_t **blen  = NULL;\n"
             "    int32_t **bexp  = NULL;\n"
             "    void **bcf      = NULL;\n"
             "    int b           = 0;\n"
             "    /* counter for randomly chosen linear forms */\n"
             "    int round = -1;\n\n"
             "    bld   = malloc(sizeof(int32_t));\n"
             "    blen  = malloc(sizeof(int32_t *));\n"
             "    bexp  = malloc(sizeof(int32_t *));\n"
             "    bcf   = malloc(sizeof(void *));\n\n"
             "    bs_t *bs    = NULL;\n"
             "    ht_t *bht   = NULL;\n"
             "    stat_t *st  = NULL;\n\n"
             "    int success = 0;\n\n"
             "    success = initialize_gba_input_data(&bs, &bht, &st,gens->lens, gens->exps, (void *)gens->cfs,	%a, 0 /* DRL order */,elim_block_len, gens->nvars,/* gens->field_char,0 [> DRL order <], gens->nvars, */ gens->ngens, saturate,	initial_hts, nr_threads, max_pairs,	update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,	info_level);\n"
             "    success = core_gba(&bs, &bht, &st);\n"
             "    if (!success) {\n"
             "      printf(\"Problem with F4, stopped computation.\\n\");\n"
             "      return 104;\n"
             "    }\n\n"
             "    export_results_from_gba(bld, blen, bexp,bcf, &malloc, &bs, &bht, &st);\n\n"
             "    int32_t *bcf_ff = (int32_t *)(*bcf);\n"
             "    int32_t *bexp_lm = get_lead_monomials(bld, blen, bexp, gens);\n"
             "    long dquot = 0;\n"
             "    int32_t *lmb= monomial_basis (bld[0], gens->nvars, bexp_lm,&dquot);\n"
             "    sp_matfglm_t  *matrix= build_matrixn(lmb, dquot, bld[0], blen, bexp, bcf_ff,bexp_lm, gens->nvars, gens->field_char);\n\n"
             "    /* display_fglm_matrix (stdout, matrix); */\n\n"
             "    if (matrix->charac != field_char) return 105;\n"
             "    if (matrix->ncols != %a) return 106;\n"
             "    if (matrix->nrows != %a) return 107;\n\n",
             str,nops(vars),char,nops(F),char,nr_col,nr_row):
    b    := is_staircase_grevlex_generic (ns,rv,lmGbxn,vars):
    cpt  := 201:
    if (b) then
        Mxn  := MultiplicationMatrix (vars[-1],ns,rv,G,ord,characteristic=char):
        print (Mxn);
        for i from 1 to nr_row do
            i2:= rv[lmGbxn[i]/vars[-1]]:
            for j from 1 to nr_col do
                fprintf (fd,
                         "    if (matrix->dense_mat[%a] != %a) return %a;\n",
                         (i-1)*nr_col+j-1,Mxn[i2,j],cpt):
                cpt:= cpt+1:
            end do:
            for j from nr_col to 1 by -1 do
                if (Mxn[i2,j] = 0) then dst[i]:= dst[i]+1:
                else break;
                end if:
            end do:
        end do:
        fprintf (fd,"\n"):
        for i from 1 to nops(triv_idx) do
            fprintf (fd,
                     "    if (matrix->triv_idx[%a] != %a) return %a;\n",
                     i-1,triv_idx[i]-1,cpt):
            cpt:= cpt+1:
        end do:
        for i from 1 to nops(triv_idx) do
            fprintf (fd,
                     "    if (matrix->triv_pos[%a] != %a) return %a;\n",
                     i-1,triv_pos[i]-1,cpt):
            cpt:= cpt+1:
        end do:
        fprintf (fd,"\n"):
        for i from 1 to nr_row do
            fprintf (fd,
                     "    if (matrix->dense_idx[%a] != %a) return %a;\n",
                     i-1,dense_idx[i]-1,cpt):
            cpt:=cpt+1:
        end do:
        fprintf (fd,"\n"):
        for i from 1 to nr_row do
            fprintf (fd,
                     "    if (matrix->dst[%a] != %a) return %a;\n",
                     i-1,dst[i],cpt):
            cpt:=cpt+1:
        end do:
        fprintf (fd,"\n"
                 "    return 0;\n}"):
    end if:
    fclose (fd):
    printf ("Writing over\n"):
end proc:
trace(main):
# F:=[x - 2*z^3 - 3*z^2 - 5*z - 7,
#     y -   z^3 +   z^2 -   z + 1,
#     z^4 -   z^3 -   z^2 -   z - 1
#    ]:
# vars:=[x,y,z]:
# char:=1073741827:
# str:="fglm_test":
# main (F,vars,char,str):
