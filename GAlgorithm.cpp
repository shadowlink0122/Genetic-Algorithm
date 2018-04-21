#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "header.h"
#include "param.h"
using namespace std;

void* my_malloc(int size){
  void* ptr = malloc(size);
  if(ptr == NULL){
    fprintf(stderr,"Faild to malloc\n");
    exit(1);
  }
  return ptr;
}

gtype_t mk_gtype(int code_length){
  gtype_t gtype = (gtype_t)my_malloc(sizeof(int) * code_length);
  return gtype;
}

void free_gtype(gtype_t gtype){
  free(gtype);
  return;
}

gtype_t mk_random_gtype(int code_length,int code_max){
  gtype_t ptr = mk_gtype(code_length);
  for(int i = 0;i < code_length;i++){
    ptr[i] = rand() % (code_max + 1);
  }
  return ptr;
}

void encode_gtype(double value,gtype_t gtype,int code_length,double minimum,double maximum){
  double gap = maximum - minimum;
  double remain_value = value - minimum;
  double value_of_code;
  int position = 1;
  int pre_code = 0;
  int i = 0;
  int tmp;
  while(i < code_length){
    value_of_code = gap / pow(2,position);
    if(remain_value >= value_of_code){
      gtype[i] = 1;
      remain_value -= value_of_code;
    }else{
      gtype[i] = 0;
    }

    if(GRAY == 1){
      tmp = gtype[i];
      gtype[i] = (pre_code) ^ (gtype[i]);
      pre_code = tmp;
    }
    position++;
    i++;
  }
  return;
}

double decode_gtype(gtype_t gtype,int code_length,double minimum,double maximum){
  double gap = maximum - minimum;
  double decoded_value = minimum;
  int position = 1;
  int pre_code = 0;
  int i = 0;
  if(GRAY == 1){
    while(i < code_length){
      pre_code = pre_code ^ gtype[i];
      if(pre_code){
        decoded_value += gap / pow(2,position);
      }
      position++;
      i++;
    }
  }else{
    while(i < code_length){
      if(gtype[i]){
        decoded_value += gap / pow(2,position);
      }
      position++;
      i++;
    }
  }
  return decoded_value;
}

void copy_gtype(gtype_t new_gtype,gtype_t old_gtype,int length){
  for(int i = 0;i < length;i++){
    new_gtype[i] = old_gtype[i];
  }
  return;
}

int cross_gtype(gtype_t gtype1,gtype_t gtype2,int length){
  int cross_point = rand() % (length - 1);
  int i = cross_point + 1;
  int tmp;
  while(i < length){
    tmp = gtype1[i];
    gtype1[i] = gtype2[i];
    gtype2[i] = tmp;
    i++;
  }
  return cross_point;
}

int mutate_gtype(gtype_t gtype,int length,int code_max,double pm){
  if(pm < 0.0 || 1.0 < pm){
    printf("%f mutation probability must be from 0.0 to 1.0 \n",pm);
    exit(-1);
  }
  int mutate_point = 0;
  double rm;
  for(int i = 0;i < length;i++){
    rm = (double)rand() / RAND_MAX;
    if(rm < pm){
      gtype[i] = rand() % (code_max + 1);
      mutate_point++;
    }
  }
  return mutate_point;
}

void print_gtype(gtype_t gtype,int length){
  printf("[");
  for(int i = 0;i < length;i++){
    if(gtype[i] < 10)printf("%d ",gtype[i]);
    else printf("(%d) ",gtype[i]);
  }
  printf("]");
}

void switch_gene(individual_t *individual){
  individual_t tmp_ptr1 = (*individual) -> next -> next;
  individual_t tmp_ptr2 = (*individual) -> next;
  (*individual) -> next -> next = *individual;
  (*individual) -> next = tmp_ptr1;
  (*individual) = tmp_ptr2;
  return;
}

individual_t mk_gene(int code_length,int code_max){
  individual_t ptr = (individual_t)malloc(sizeof(ga_individual));
  ptr->gtype = mk_random_gtype(code_length,code_max);
  ptr->ptype = 0;
  ptr->fitness = 0;
  ptr->next = NULL;
  ptr->parent1 = 0;
  ptr->parent2 = 0;
  ptr->cross_point = 0;
  return ptr;
}

void copy_gene(individual_t new_gene,individual_t old_gene,int code_length){
  copy_gtype(new_gene->gtype,old_gene->gtype,code_length);
  new_gene->ptype = old_gene->ptype;
  new_gene->fitness = old_gene->fitness;
  new_gene->parent1 = old_gene->rank;
  new_gene->parent2 = old_gene->rank;
  new_gene->cross_point = code_length - 1;
  return;
}

int mk_children_genes(individual_t child1,individual_t child2,individual_t parent1,individual_t parent2,int code_length,int code_max,double pm){
  int cross_point,mutateCount;
  copy_gene(child1,parent1,code_length);
  copy_gene(child2,parent2,code_length);
  cross_point = cross_gtype(child1->gtype,child2->gtype,code_length);
  child1->parent1 = parent1->rank;
  child1->parent2 = parent2->rank;
  child1->cross_point = cross_point;
  child2->parent1 = parent2->rank;
  child2->parent2 = parent1->rank;
  child2->cross_point = cross_point;
  mutateCount = mutate_gtype(child1->gtype,code_length,code_max,pm);
  mutateCount += mutate_gtype(child2->gtype,code_length,code_max,pm);
  return mutateCount;
}

ga_population_t mk_init_ga_population(int population_size,int code_length,int code_max){
  ga_population_t population = (ga_population_t)malloc(sizeof(ga_population));
  population->pselect = (double*)my_malloc(sizeof(double) * population_size);
  population->mutate_count = 0;
  population->population_size = population_size;
  population->code_length = code_length;
  population->code_max = code_max;
  individual_t list_tale;
  population->genes = mk_gene(code_length,code_max);
  list_tale = population->genes;
  for(int i = 1;i < population_size;i++){
    list_tale->next = mk_gene(code_length,code_max);
    list_tale = list_tale->next;
  }
  return population;
}

void print_sequence(char ch,int length){
  for(int i = 0;i < length;i++)printf("%c",ch);
}

void print_population(ga_population_t population){
  individual_t member = population->genes;
  int i = 0;
  printf("---------");
  print_sequence('-',LENGTH + 2);
  printf("-------\n");
  printf("# parents xsite gtype");
  print_sequence('-',LENGTH - 3);
  printf("ptype fitness\n");

  while(member != NULL){
    printf("%-3d (%3d,%3d) %3d ",i++,member->parent1,member->parent2,member->cross_point);
    print_gtype(member->gtype,population->code_length);
    printf(" %+3.3f %+3.3f\n",member->ptype,member->fitness);
    member = member->next;
  }
  printf("total mutate %d\n",population->mutate_count);
  return;
}

void print_fitness(ga_population_t population){
  printf("%f, %f, %f, %f, ",population->max_fitness,population->avg_fitness,population->min_fitness,population->genes->ptype);
  print_gtype(population->genes->gtype,population->code_length);
  printf("\n");
  return;
}

int less_than(individual_t individualA,individual_t individualB){
  return (individualA->fitness < individualB->fitness);
}

void calc_fitness(ga_population_t population,double value_min,double value_max){
  individual_t ptr = population->genes;
  individual_t next;
  individual_t individual_ptr = NULL;
  individual_t search_ptr = ptr;
  double x,y;
  while(ptr != NULL){
    x = decode_gtype(ptr->gtype,population->code_length,value_min,value_max);
    ptr->ptype = x;
    y = F_X(x);
    ptr->fitness = G_Y(y);
    next = ptr->next;
    search_ptr = individual_ptr;
    if(search_ptr == NULL || less_than(individual_ptr,ptr)){
      ptr->next = individual_ptr;
      individual_ptr = ptr;
    }else{
      while(search_ptr->next != NULL){
        if(less_than(search_ptr->next,ptr))break;
        search_ptr = search_ptr->next;
      }
      ptr->next = search_ptr->next;
      search_ptr->next = ptr;
    }
    ptr = next;
  }
  population->genes = individual_ptr;
  return;
}

void calc_pselect(ga_population_t population){
  population->pselect[0] = population->genes->fitness;
  individual_t genes_ptr = population->genes->next;
  for(int i = 1;i < population->population_size;i++){
    population->pselect[i] = population->pselect[i-1] + genes_ptr->fitness;
    genes_ptr = genes_ptr->next;
  }
  for(int i = 0;i < population->population_size;i++){
    population->pselect[i] /= population->pselect[population->population_size-1];
  }
}

individual_t select_parent_roulette(ga_population_t population){
  int i = 0;
  double r;
  individual_t parent;
  r = (double)rand() / RAND_MAX;
  parent = population->genes;
  while(r > population->pselect[i++])parent = parent->next;
  return parent;
}

individual_t select_parent_tournament(ga_population_t population,int tournament_size){
  int pop_t = population->population_size;
  int r,min_t = pop_t;
  individual_t min_selected = NULL;
  individual_t ptr;
  for(int i = 0;i < tournament_size;i++){
    r = rand() % pop_t;
    if(min_t > r)min_t = r;
  }
  ptr = population->genes;
  for(int i = 0;i < min_t;i++)ptr = ptr->next;
  min_selected = ptr;
  return min_selected;
}

individual_t select_parent(ga_population_t population){
  individual_t parent;
  switch(SELECTION_METHOD){
  case 1:
    parent = select_parent_roulette(population);
    break;
  case 2:
    parent = select_parent_tournament(population,TOURNAMENT_SIZE);
    break;
  default:
    fprintf(stderr,"invalid member on SELECTION_METHOD\n");
    exit(1);
  }
  return parent;
}

void normalize_population(ga_population_t population){
  individual_t tmp = population->genes;
  population->max_fitness = population->genes->fitness;
  double avg = 0.0;
  for(int i = 0;i < population->population_size;i++){
    avg += tmp->fitness;
    tmp->rank = i;
    if(tmp->next == NULL){
      population->min_fitness = tmp->fitness;
    }
    tmp = tmp->next;  
  }
  avg = avg / population->population_size;
  population->avg_fitness = avg;
  return;
}

void generate_population(ga_population_t new_population,ga_population_t old_population,double gap,double elite_rate,double mutate_prob,double crossover_prob){
  int num_of_remain = (int)(old_population->population_size * (1 - gap));
  int num_of_elite = (num_of_remain * elite_rate);
  int generated;
  double rand_double;
  individual_t old_gene = old_population->genes;
  individual_t new_gene = new_population->genes;
  calc_pselect(old_population);
  for(generated = 0;generated < num_of_elite;generated++){
    copy_gene(new_gene,old_gene,old_population->code_length);
    old_gene = old_gene->next;
    new_gene = new_gene->next;
  }
  for(;generated < num_of_remain;generated++){
    copy_gene(new_gene,select_parent(old_population),old_population->code_length);
    new_gene = new_gene->next;
  }
  new_population->mutate_count = 0;
  if((old_population->population_size - generated) % 2){
    copy_gene(new_gene,select_parent(old_population),old_population->code_length);
    new_population->mutate_count += mutate_gtype(new_gene->gtype,old_population->code_length,old_population->code_max,mutate_prob);
    new_gene = new_gene->next;
    generated++;
  }
  for(;generated < old_population->population_size;generated += 2){
    rand_double = (double)rand() / RAND_MAX;
    if(rand_double < crossover_prob){
      new_population->mutate_count += mk_children_genes(new_gene,new_gene->next,select_parent(old_population),select_parent(old_population),old_population->code_length,old_population->code_max,mutate_prob);
      new_gene = new_gene->next->next;
    }else{
      copy_gene(new_gene,select_parent(old_population),old_population->code_length);
      new_population->mutate_count += mutate_gtype(new_gene->gtype,old_population->code_length,old_population->code_max,mutate_prob);
      new_gene = new_gene->next;
      copy_gene(new_gene,select_parent(old_population),old_population->code_length);
      new_population->mutate_count += mutate_gtype(new_gene->gtype,old_population->code_length,old_population->code_max,mutate_prob);
      new_gene = new_gene->next;
    }
  }
  return;
}

int main_ga(){
  srand(time(NULL));
  ga_population_t parent_group = mk_init_ga_population(POP,LENGTH,CODE_MAX);
  ga_population_t child_group = mk_init_ga_population(POP,LENGTH,CODE_MAX);
  if(PRINT_FITNESS == 1)
    printf("#generation,max_fitness,avg_fitness,min_fitness,beet_individual_ptype,beet_individual_ptype\n");
  
  for(int i = 0;i <= GENERATION;i++){
    calc_fitness(parent_group,MIN,MAX);
    normalize_population(parent_group);
    if(PRINT_GROUP == 1)print_population(parent_group);
    if(PRINT_FITNESS == 1){
      printf("%3d, ",i);
      print_fitness(parent_group);
    }
    generate_population(child_group,parent_group,GAP,ELITE_RATE,P_MUTATE,P_CROSS);
    parent_group = child_group;
  }
  return 0;
}

int main(){
  main_ga();
  return 0;
}