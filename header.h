
typedef int* gtype_t;
typedef struct ga_individual* individual_t;
typedef struct ga_population* ga_population_t;

struct ga_individual{
  gtype_t gtype;
  double ptype;
  double fitness;
  individual_t next;
  int rank;
  int parent1;
  int parent2;
  int cross_point;
};

struct ga_population{
  individual_t genes;
  double *pselect;
  int mutate_count;
  double max_fitness;
  double min_fitness;
  double avg_fitness;
  int population_size;
  int code_length;
  int code_max;
};